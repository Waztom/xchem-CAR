import React, { useCallback, useLayoutEffect } from 'react';
import {
  colors,
  Divider,
  makeStyles,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  TableSortLabel,
  Tooltip
} from '@material-ui/core';
import { useTable, useSortBy, useExpanded, useRowSelect, useFilters } from 'react-table';
import { TargetRow } from './components/TargetRow';
import { setFilters, useBatchesTableStateStore } from '../../../common/stores/batchesTableStateStore';
import { useBatchContext } from '../../hooks/useBatchContext';
import { TableToolbar } from '../TableToolbar';
import { useTableSelectionColumn } from './hooks/useTableSelectionColumn';
import { useTableColumns } from './hooks/useTableColumns';
import { useGetTableData } from '../../hooks/useGetTableData';

const useStyles = makeStyles(theme => ({
  root: {
    display: 'grid',
    width: '100%'
  },
  table: {
    display: 'grid',
    overflowX: 'auto',
    '& tr': {
      display: 'grid',
      alignItems: 'stretch',
      gap: `0 ${theme.spacing()}px`,
      borderBottom: `1px solid ${theme.palette.divider}`,
      paddingLeft: theme.spacing(2),
      paddingRight: theme.spacing(2)
    },
    '& th, td': {
      display: 'grid',
      placeItems: 'center',
      border: 0,
      padding: 0
    }
  },
  flexCell: {
    display: 'flex',
    justifyContent: 'center',
    alignItems: 'center',
    gap: theme.spacing(1 / 2)
  },
  row: {
    gridTemplateColumns: ({ maxNoSteps }) => {
      if (!maxNoSteps) {
        return '40px 60px';
      }
      return `40px 60px repeat(${maxNoSteps}, 282px)`;
    }
  },
  sortIconInactive: {
    '& svg': {
      color: `${colors.grey[400]} !important`
    }
  }
}));

const getRowId = (row, relativeIndex, parent) => {
  return parent ? [parent.id, row.id].join('.') : String(row.id);
};

export const BatchTable = () => {
  const tableData = useGetTableData();

  const batch = useBatchContext();

  const expanded = useBatchesTableStateStore(useCallback(state => state.expanded[batch.id] || {}, [batch.id]));
  const selected = useBatchesTableStateStore(useCallback(state => state.selected[batch.id] || {}, [batch.id]));
  const filters = useBatchesTableStateStore(useCallback(state => state.filters[batch.id] || [], [batch.id]));

  const maxNoSteps = tableData.length
    ? Math.max(
        ...tableData
          .map(({ subRows }) => {
            if (subRows.length) {
              return subRows.map(({ reactions }) => reactions.length);
            }
            return 0;
          })
          .flat()
      )
    : 0;

  const columns = useTableColumns(maxNoSteps);

  const tableInstance = useTable(
    {
      columns,
      data: tableData,
      getRowId,
      initialState: {
        expanded,
        selectedRowIds: selected,
        filters,
        paginateExpandedRows: true,
        hiddenColumns: [
          'catalogentries',
          'otchem',
          'target-vendor',
          'target-leadtime',
          'target-price',
          ...new Array(maxNoSteps)
            .fill(0)
            .map((_, index) => [
              `reactant-vendor-step-${index}`,
              `reactant-preferred-vendor-step-${index}`,
              `reactant-preferred-leadtime-step-${index}`,
              `reactant-preferred-price-step-${index}`
            ])
            .flat(),
          'reactant-smiles'
        ]
      }
    },
    useFilters,
    useSortBy,
    useExpanded,
    useRowSelect,
    useTableSelectionColumn
  );

  const { getTableProps, getTableBodyProps, headerGroups, prepareRow, rows, state } = tableInstance;

  useLayoutEffect(() => {
    // Sync filters with zustand store
    setFilters(batch.id, state.filters);
  }, [batch.id, state.filters]);

  const classes = useStyles({ maxNoSteps });

  return (
    <div className={classes.root}>
      <TableToolbar tableInstance={tableInstance} />
      <Divider />
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
          {rows
            // When filters are applied, there might be a case where children partially match multiple filters but not all.
            // In that case react-library still displays the parent row even though there are no subRows.
            .filter(row =>
              row.depth === 0 && row.subRows.length !== row.originalSubRows.length ? row.subRows.length : true
            )
            .map(row => {
              prepareRow(row);

              if (row.depth === 0) {
                return <TargetRow {...row.getRowProps()} row={row} />;
              }

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
    </div>
  );
};
