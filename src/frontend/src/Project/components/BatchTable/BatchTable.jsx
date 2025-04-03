import React, { useCallback, useLayoutEffect } from 'react';
import {
  colors,
  Divider,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  TableSortLabel,
  Tooltip,
} from '@mui/material';
import { styled } from '@mui/material/styles';
import { useTable, useSortBy, useExpanded, useRowSelect, useFilters } from 'react-table';
import { SuspenseWithBoundary } from '../../../common/components/SuspenseWithBoundary';
import { TargetRow } from './components/TargetRow';
import { setFilters, useBatchesTableStateStore } from '../../../common/stores/batchesTableStateStore';
import { useBatchContext } from '../../hooks/useBatchContext';
import { TableToolbar } from '../TableToolbar';
import { useTableSelectionColumn } from './hooks/useTableSelectionColumn';
import { useTableColumns } from './hooks/useTableColumns';
import { useGetTableData } from '../../hooks/useGetTableData';

const TableContainer = styled('div')(({ theme }) => ({
  display: 'grid',
  width: '100%'
}));

const StyledTable = styled(Table)(({ theme }) => ({
  display: 'grid',
  overflowX: 'auto',
  '& tr': {
    display: 'grid',
    alignItems: 'stretch',
    gap: theme.spacing(1),
    borderBottom: `1px solid ${theme.palette.divider}`,
    padding: theme.spacing(0, 2)
  },
  '& th, td': {
    display: 'grid',
    placeItems: 'center',
    border: 0,
    padding: 0
  }
}));

const FlexCell = styled('div')(({ theme }) => ({
  display: 'flex',
  justifyContent: 'center',
  alignItems: 'center',
  gap: theme.spacing(1/2)
}));

const StyledTableRow = styled(TableRow, {
  shouldForwardProp: prop => prop !== 'maxNoSteps'
})(({ theme, maxNoSteps }) => ({
  gridTemplateColumns: maxNoSteps 
    ? `40px 60px repeat(${maxNoSteps}, 282px)`
    : '40px 60px'
}));

const StyledSortLabel = styled(TableSortLabel)(({ theme, active }) => ({
  '& svg': {
    color: !active ? colors.grey[400] : undefined
  }
}));

const BatchTableContent = () => {
  const tableData = useGetTableData();
  const batch = useBatchContext();
  const expanded = useBatchesTableStateStore(useCallback(state => state.expanded[batch.id] || {}, [batch.id]));
  const selected = useBatchesTableStateStore(useCallback(state => state.selected[batch.id] || {}, [batch.id]));
  const filters = useBatchesTableStateStore(useCallback(state => state.filters[batch.id] || [], [batch.id]));

  const maxNoSteps = tableData.length
    ? Math.max(...tableData.map(({ subRows }) => 
        subRows.length ? subRows.map(({ reactions }) => reactions.length) : 0).flat())
    : 0;

  const columns = useTableColumns(maxNoSteps);
  const tableInstance = useTable(
    {
      columns,
      data: tableData,
      getRowId: (row, relativeIndex, parent) => parent ? [parent.id, row.id].join('.') : String(row.id),
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
    setFilters(batch.id, state.filters);
  }, [batch.id, state.filters]);

  return (
    <TableContainer>
      <SuspenseWithBoundary>
        <TableToolbar tableInstance={tableInstance} />
      </SuspenseWithBoundary>
      <Divider />
      <StyledTable {...getTableProps()}>
        <TableHead>
          {headerGroups.map(headerGroup => {
            const { key: headerGroupKey, ...headerGroupProps } = headerGroup.getHeaderGroupProps();
            return (
              <StyledTableRow 
                key={headerGroupKey}
                {...headerGroupProps} 
                maxNoSteps={maxNoSteps}
              >
                {headerGroup.headers.map(column => {
                  const { key: headerKey, ...headerProps } = column.getHeaderProps();
                  if (column.canSort) {
                    const { title, ...sortProps } = column.getSortByToggleProps();
                    return (
                      <Tooltip key={headerKey} title={`Sort by ${column.sortLabel}`} {...headerProps}>
                        <TableCell {...sortProps}>
                          <FlexCell>
                            {column.render('Header')}
                            <StyledSortLabel
                              active={true}
                              direction={column.isSortedDesc ? 'desc' : 'asc'}
                            />
                          </FlexCell>
                        </TableCell>
                      </Tooltip>
                    );
                  }
                  return (
                    <TableCell key={headerKey} {...headerProps}>
                      {column.render('Header')}
                    </TableCell>
                  );
                })}
              </StyledTableRow>
            );
          })}
        </TableHead>
        <TableBody {...getTableBodyProps()}>
          {rows
            .filter(row => row.depth === 0 && row.subRows.length !== row.originalSubRows.length 
              ? row.subRows.length 
              : true
            )
            .map(row => {
              prepareRow(row);
              const { key: rowKey, ...rowProps } = row.getRowProps();
              
              if (row.depth === 0) {
                return (
                  <SuspenseWithBoundary key={rowKey}>
                    <TargetRow {...rowProps} row={row} />
                  </SuspenseWithBoundary>
                );
              }
              
              return (
                <StyledTableRow 
                  key={rowKey}
                  {...rowProps} 
                  maxNoSteps={maxNoSteps}
                >
                  {row.cells.map(cell => {
                    const { key: cellKey, ...cellProps } = cell.getCellProps();
                    return (
                      <TableCell key={cellKey} {...cellProps}>
                        {cell.render('Cell')}
                      </TableCell>
                    );
                  })}
                </StyledTableRow>
              );
            })}
        </TableBody>
      </StyledTable>
    </TableContainer>
  );
};

export const BatchTable = () => (
  <SuspenseWithBoundary>
    <BatchTableContent />
  </SuspenseWithBoundary>
);

BatchTable.displayName = 'BatchTable';
