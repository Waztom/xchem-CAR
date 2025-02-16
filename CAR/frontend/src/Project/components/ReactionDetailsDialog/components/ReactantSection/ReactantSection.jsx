import React from 'react';
import {
  colors,
  makeStyles,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  TableSortLabel,
  Tooltip,
  Typography
} from '@material-ui/core';
import { DialogSection } from '../../../../../common/components/DialogSection';
import { DialogSectionHeading } from '../../../../../common/components/DialogSectionHeading';
import { useReactantProductTableColumns } from './hooks/useReactantProductTableColumns';
import { useTable, useSortBy } from 'react-table';
import PubChem from '../../../../../assets/pubchem.svg';
import PubChemSafety from '../../../../../assets/pubchem-safety.svg';

const useStyles = makeStyles(theme => ({
  table: {
    '& tr': {
      display: 'grid',
      alignItems: 'stretch',
      gap: `0 ${theme.spacing()}px`,
      borderBottom: `1px solid ${theme.palette.divider}`,
      paddingLeft: theme.spacing(2),
      paddingRight: theme.spacing(2),
      gridTemplateColumns: '160px 1fr 105px repeat(2, 95px)'
    },
    '& th, td': {
      display: 'grid',
      placeItems: 'center start',
      border: 0,
      padding: 0
    }
  },
  flexCell: {
    display: 'flex',
    alignItems: 'center',
    gap: theme.spacing(1 / 2)
  },
  sortIconInactive: {
    '& svg': {
      color: `${colors.grey[400]} !important`
    }
  },
  sections: {
    display: 'grid',
    gridTemplateColumns: 'repeat(3, 1fr)'
  }
}));

export const ReactantSection = ({ reactant, index }) => {
  const classes = useStyles();

  const reactantpubcheminfo = reactant.reactantpubcheminfo || {};

  const columns = useReactantProductTableColumns();

  const { getTableProps, headerGroups, getTableBodyProps, rows, prepareRow } = useTable(
    {
      columns,
      data: reactant.catalogentries
    },
    useSortBy
  );

  return (
    <DialogSection>
      <DialogSectionHeading>Reactant {index + 1}</DialogSectionHeading>
      <div className={classes.sections}>
        <div>
          <Typography>
            Smiles: <strong>{reactant.smiles}</strong>
          </Typography>
          {reactantpubcheminfo.cas && (
            <Typography>
              CAS number: <strong>{reactantpubcheminfo.cas}</strong>
            </Typography>
          )}
        </div>

        {!!reactantpubcheminfo.summaryurl && (
          <div>
            <Tooltip title="PubChem compound summary">
              <a href={reactantpubcheminfo.summaryurl} target="_blank" rel="noreferrer">
                <img src={PubChem} height={50} />
              </a>
            </Tooltip>
          </div>
        )}
        {!!reactantpubcheminfo.lcssurl && (
          <div>
            <Tooltip title="Laboratory chemical safety summary">
              <a href={reactantpubcheminfo.lcssurl} target="_blank" rel="noreferrer">
                <img src={PubChemSafety} height={50} />
              </a>
            </Tooltip>
          </div>
        )}
      </div>

      {!!reactant.catalogentries.length && (
        <>
          <Typography>Catalog information: </Typography>
          <Table className={classes.table} {...getTableProps()}>
            <TableHead>
              {headerGroups.map(headerGroup => (
                <TableRow {...headerGroup.getHeaderGroupProps()} className={classes.row}>
                  {headerGroup.headers.map(column => {
                    if (column.canSort) {
                      // Title is unused
                      const { title, ...rest } = column.getSortByToggleProps();

                      return (
                        <Tooltip title={`Sort by ${column.sortLabel}`} {...column.getHeaderProps()}>
                          <TableCell {...rest}>
                            <div className={classes.flexCell}>
                              {column.render('Header')}
                              <TableSortLabel
                                className={!column.isSorted ? classes.sortIconInactive : undefined}
                                active={true}
                                direction={column.isSortedDesc ? 'desc' : 'asc'}
                              />
                            </div>
                          </TableCell>
                        </Tooltip>
                      );
                    }
                    return <TableCell {...column.getHeaderProps()}>{column.render('Header')}</TableCell>;
                  })}
                </TableRow>
              ))}
            </TableHead>
            <TableBody {...getTableBodyProps()}>
              {rows.map(row => {
                prepareRow(row);

                return (
                  <TableRow {...row.getRowProps()} className={classes.row}>
                    {row.cells.map(cell => (
                      <TableCell {...cell.getCellProps()}>{cell.render('Cell')}</TableCell>
                    ))}
                  </TableRow>
                );
              })}
            </TableBody>
          </Table>
        </>
      )}
    </DialogSection>
  );
};
