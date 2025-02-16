import React, { useMemo } from 'react';
import { colors, makeStyles, Tooltip, Typography } from '@material-ui/core';
import { IconComponent } from '../../../../common/components/IconComponent';
import { IoFootsteps } from 'react-icons/io5';
import { AutocompleteFilter } from '../components/AutocompleteFilter';
import {
  createTableMethodAutocompleteFilter,
  createTableMethodYesNoFilter,
  createTableMethodSmilesFilter
} from '../utils/createTableMethodFilters';
import {
  createTableTargetAutocompleteFilter,
  createTableTargetRangeFilter,
  createTableTargetYesNoFilter
} from '../utils/createTableTargetFilters';
import { YesNoFilter } from '../components/YesNoFilter';
import { RangeFilter } from '../components/RangeFilter';
import { PREFERRED_LEAD_TIME, PREFERRED_PRICE, PREFERRED_VENDORS } from '../../../constants/preferredContstants';
import { PreferredFlagIndicator } from '../components/PreferredFlagIndicator/PreferredFlagIndicator';
import { formatPreferredVendorsString } from '../../../utils/formatPreferredVendorsString';
import { SmilesFilter } from '../components/SmilesFilter/SmilesFilter';
import classNames from 'classnames';
import { FormatListNumbered } from '@material-ui/icons';
import { requestReactionDetailsDialog } from '../../../stores/reactionDetailsDialogStore';

const useStyles = makeStyles(theme => ({
  text: {
    width: '100%',
    textAlign: 'center'
  },
  flexCell: {
    display: 'flex',
    justifyContent: 'center',
    alignItems: 'center',
    gap: theme.spacing(1 / 2)
  },
  reactionWrapper: {
    display: 'grid',
    cursor: 'pointer'
  },
  reactionNameWrapper: {
    textAlign: 'center',
    display: 'grid',
    lineHeight: 1
  },
  preferredIndicatorsWrapper: {
    display: 'grid',
    gap: theme.spacing()
  },
  unsuccessful: {
    backgroundColor: colors.red[100]
  },
  image: {
    mixBlendMode: 'multiply'
  }
}));

const filterByMethodReactionName = createTableMethodAutocompleteFilter((row, ids, filterValue) =>
  filterValue.includes(row.values[ids[0]])
);

const filterByMethodOTExecutable = createTableMethodYesNoFilter(
  (row, ids, filterValue) => row.values[ids[0]] === filterValue
);

const filterByTargetCatalogEntry = createTableTargetYesNoFilter(
  (row, ids, filterValue) => !!row.values[ids[0]].length === filterValue
);

const filterByTargetVendor = createTableTargetAutocompleteFilter((row, ids, filterValue) => {
  const vendors = row.original.catalogentries.map(({ vendor }) => vendor);
  return filterValue.some(value => vendors.includes(value));
});

const filterByTargetLeadTime = createTableTargetRangeFilter((row, ids, filterValue) => {
  const leadTimes = row.original.catalogentries.map(({ leadtime }) => leadtime);
  const [min, max] = filterValue;
  return leadTimes.some(leadTime => min >= leadTime && leadTime <= max);
});

const filterByTargetPrice = createTableTargetRangeFilter((row, ids, filterValue) => {
  const prices = row.original.catalogentries.map(({ upperprice }) => upperprice);
  const [min, max] = filterValue;
  return prices.some(price => price >= min && price <= max);
});

const filterByMethodReactantVendor = index =>
  createTableMethodAutocompleteFilter((row, ids, filterValue) => {
    const vendors = (
      row.original.reactions?.[index]?.reactants?.map(({ catalogentries }) =>
        catalogentries?.map(({ vendor }) => vendor)
      ) || []
    ).flat(2);
    return filterValue.some(value => vendors.includes(value));
  });

const filterByMethodPreferredReactantVendor = index =>
  createTableMethodYesNoFilter(
    (row, ids, filterValue) => filterValue === row.original.reactions?.[index]?.preferredVendor
  );

const filterByMethodPreferredReactantLeadTime = index =>
  createTableMethodYesNoFilter(
    (row, ids, filterValue) => filterValue === row.original.reactions?.[index]?.preferredLeadTime
  );

const filterByMethodPreferredReactantPrice = index =>
  createTableMethodYesNoFilter(
    (row, ids, filterValue) => filterValue === row.original.reactions?.[index]?.preferredPrice
  );

const filterByMethodReactantsExcludeSmiles = createTableMethodSmilesFilter((row, ids, filterValue) => {
  const smiles = (row.original.reactions?.map(({ reactants }) => reactants?.map(({ smiles }) => smiles)) || []).flat();
  return !smiles.some(smile => filterValue.includes(smile));
});

export const useTableColumns = maxNoSteps => {
  const classes = useStyles();

  // Since filter functions should be memoized, we can't call them directly
  const reactantVendorFilters = useMemo(() => {
    return new Array(maxNoSteps).fill(0).map((_, index) => filterByMethodReactantVendor(index));
  }, [maxNoSteps]);
  const preferredReactantVendorFilters = useMemo(() => {
    return new Array(maxNoSteps).fill(0).map((_, index) => filterByMethodPreferredReactantVendor(index));
  }, [maxNoSteps]);
  const preferredReactantLeadTimeFilters = useMemo(() => {
    return new Array(maxNoSteps).fill(0).map((_, index) => filterByMethodPreferredReactantLeadTime(index));
  }, [maxNoSteps]);
  const preferredReactantPriceFilters = useMemo(() => {
    return new Array(maxNoSteps).fill(0).map((_, index) => filterByMethodPreferredReactantPrice(index));
  }, [maxNoSteps]);

  const columns = useMemo(() => {
    return [
      {
        accessor: 'position',
        sortLabel: 'position',
        disableFilters: true,
        Header: () => {
          return <FormatListNumbered />;
        },
        Cell: ({ value }) => {
          return (
            <Typography className={classes.text} component="span" noWrap>
              {value}
            </Typography>
          );
        },
        sortType: 'number'
      },
      ...new Array(maxNoSteps).fill(0).map((_, index) => {
        return {
          accessor: `reactions[${index}].reactionclass`,
          sortLabel: `reaction step ${index + 1}`,
          filterOrder: 6,
          Header: () => {
            return (
              <div className={classes.flexCell}>
                <IconComponent Component={IoFootsteps} />
                <Typography component="span">{index + 1}</Typography>
              </div>
            );
          },
          Cell: ({ value, row }) => {
            const reaction = row.original.reactions[index];

            if (!reaction) {
              return null;
            }

            return (
              <div className={classNames(classes.flexCell, !reaction.success && classes.unsuccessful)}>
                <Tooltip title="Show reaction details">
                  <div className={classes.reactionWrapper} onClick={() => requestReactionDetailsDialog(reaction)}>
                    <img className={classes.image} src={reaction.image} width={270} height={60} />
                    <div className={classes.reactionNameWrapper}>
                      <Typography variant="caption" noWrap>
                        {value}
                      </Typography>
                    </div>
                  </div>
                </Tooltip>
                <div className={classes.preferredIndicatorsWrapper}>
                  <PreferredFlagIndicator reaction={reaction} type="vendor" />
                  <PreferredFlagIndicator reaction={reaction} type="leadTime" />
                  <PreferredFlagIndicator reaction={reaction} type="price" />
                </div>
              </div>
            );
          },
          Filter: ({ column: { filterValue, setFilter }, preFilteredFlatRows }) => {
            return (
              <AutocompleteFilter
                id={`reaction-${index + 1}-filter`}
                options={[
                  ...new Set(
                    preFilteredFlatRows
                      .filter(row => row.depth === 1)
                      .map(row => row.original.reactions?.[index]?.reactionclass)
                  )
                ].sort()}
                label={`Reaction type - step ${index + 1}`}
                placeholder="Reaction type"
                filterValue={filterValue}
                setFilter={setFilter}
              />
            );
          },
          filter: filterByMethodReactionName
        };
      }),
      {
        accessor: 'otchem',
        filterOrder: 1,
        Filter: ({ column: { filterValue, setFilter } }) => {
          return (
            <YesNoFilter
              id="otchem-filter"
              label="Executable on OpenTrons"
              filterValue={filterValue}
              setFilter={setFilter}
            />
          );
        },
        filter: filterByMethodOTExecutable
      },
      {
        accessor: 'catalogentries',
        filterOrder: 2,
        Filter: ({ column: { filterValue, setFilter } }) => {
          return (
            <YesNoFilter
              id="target-catalogentries-filter"
              label="Target has catalog entry"
              filterValue={filterValue}
              setFilter={setFilter}
            />
          );
        },
        filter: filterByTargetCatalogEntry
      },
      {
        id: 'target-vendor',
        defaultCanFilter: true,
        filterOrder: 3,
        Filter: ({ column: { filterValue, setFilter }, preFilteredFlatRows }) => {
          return (
            <AutocompleteFilter
              id="target-vendor-filter"
              options={[
                ...new Set(
                  preFilteredFlatRows
                    .filter(row => row.depth === 0)
                    .map(row => row.original.catalogentries?.map(({ vendor }) => vendor) || [])
                    .flat()
                )
              ].sort()}
              label="Target vendor"
              placeholder="Target vendor"
              filterValue={filterValue}
              setFilter={setFilter}
            />
          );
        },
        filter: filterByTargetVendor
      },
      {
        id: 'target-leadtime',
        defaultCanFilter: true,
        filterOrder: 4,
        Filter: ({ column: { filterValue, setFilter }, preFilteredFlatRows }) => {
          const leadTimes = preFilteredFlatRows
            .filter(row => row.depth === 0)
            .map(row => row.original.catalogentries?.map(({ leadtime }) => leadtime) || [])
            .flat();
          return (
            <RangeFilter
              id="target-leadtime"
              label="Target lead time"
              min={Math.min(...leadTimes)}
              max={Math.max(...leadTimes)}
              filterValue={filterValue}
              setFilter={setFilter}
            />
          );
        },
        filter: filterByTargetLeadTime
      },
      {
        id: 'target-price',
        defaultCanFilter: true,
        filterOrder: 5,
        Filter: ({ column: { filterValue, setFilter }, preFilteredFlatRows }) => {
          const prices = preFilteredFlatRows
            .filter(row => row.depth === 0)
            .map(row => row.original.catalogentries?.map(({ upperprice }) => upperprice) || [])
            .flat();

          return (
            <RangeFilter
              id="target-price"
              label="Target price"
              min={Math.min(...prices)}
              max={Math.max(...prices)}
              filterValue={filterValue}
              setFilter={setFilter}
            />
          );
        },
        filter: filterByTargetPrice
      },
      ...new Array(maxNoSteps).fill(0).map((_, index) => {
        return {
          id: `reactant-vendor-step-${index}`,
          defaultCanFilter: true,
          filterOrder: 7,
          Filter: ({ column: { filterValue, setFilter }, preFilteredFlatRows }) => {
            return (
              <AutocompleteFilter
                id={`reactant-vendor-${index + 1}-filter`}
                options={[
                  ...new Set(
                    preFilteredFlatRows
                      .filter(row => row.depth === 1)
                      .map(row =>
                        row.original.reactions?.[index]?.reactants.map(({ catalogentries }) =>
                          catalogentries.map(({ vendor }) => vendor)
                        )
                      )
                      .flat(2)
                      .filter(vendor => !!vendor)
                  )
                ].sort()}
                label={`Reactant vendor - step ${index + 1}`}
                placeholder="Reactant vendor"
                filterValue={filterValue}
                setFilter={setFilter}
              />
            );
          },
          filter: reactantVendorFilters[index]
        };
      }),
      ...new Array(maxNoSteps).fill(0).map((_, index) => {
        return {
          id: `reactant-preferred-vendor-step-${index}`,
          defaultCanFilter: true,
          filterOrder: 8,
          Filter: ({ column: { filterValue, setFilter } }) => {
            return (
              <YesNoFilter
                id={`reactant-preferred-vendor-${index + 1}-filter`}
                label={`Reactant vendor is ${formatPreferredVendorsString(PREFERRED_VENDORS)} - step ${index + 1}`}
                filterValue={filterValue}
                setFilter={setFilter}
              />
            );
          },
          filter: preferredReactantVendorFilters[index]
        };
      }),
      ...new Array(maxNoSteps).fill(0).map((_, index) => {
        return {
          id: `reactant-preferred-leadtime-step-${index}`,
          defaultCanFilter: true,
          filterOrder: 9,
          Filter: ({ column: { filterValue, setFilter } }) => {
            return (
              <YesNoFilter
                id={`reactant-preferred-leadtime-${index + 1}-filter`}
                label={`Reactant lead time within ${PREFERRED_LEAD_TIME} weeks - step ${index + 1}`}
                filterValue={filterValue}
                setFilter={setFilter}
              />
            );
          },
          filter: preferredReactantLeadTimeFilters[index]
        };
      }),
      ...new Array(maxNoSteps).fill(0).map((_, index) => {
        return {
          id: `reactant-preferred-price-step-${index}`,
          defaultCanFilter: true,
          filterOrder: 10,
          Filter: ({ column: { filterValue, setFilter } }) => {
            return (
              <YesNoFilter
                id={`reactant-preferred-price-${index + 1}-filter`}
                label={`Reactant price within ${PREFERRED_PRICE} - step ${index + 1}`}
                filterValue={filterValue}
                setFilter={setFilter}
              />
            );
          },
          filter: preferredReactantPriceFilters[index]
        };
      }),
      {
        id: 'reactant-smiles',
        defaultCanFilter: true,
        filterOrder: 11,
        Filter: ({ column: { filterValue, setFilter } }) => {
          return (
            <SmilesFilter
              id="reactant-smiles"
              label="Exclude reactant smiles"
              filterValue={filterValue}
              setFilter={setFilter}
            />
          );
        },
        filter: filterByMethodReactantsExcludeSmiles
      }
    ];
  }, [
    maxNoSteps,
    classes.text,
    classes.flexCell,
    classes.unsuccessful,
    classes.reactionWrapper,
    classes.image,
    classes.reactionNameWrapper,
    classes.preferredIndicatorsWrapper,
    reactantVendorFilters,
    preferredReactantVendorFilters,
    preferredReactantLeadTimeFilters,
    preferredReactantPriceFilters
  ]);

  return columns;
};
