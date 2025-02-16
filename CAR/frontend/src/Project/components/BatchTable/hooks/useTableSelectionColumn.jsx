import React from 'react';
import { Checkbox } from '@material-ui/core';
import { useBatchContext } from '../../../hooks/useBatchContext';
import { setRowsSelected } from '../../../../common/stores/batchesTableStateStore';

export const useTableSelectionColumn = hooks => {
  const batch = useBatchContext();

  hooks.visibleColumns.push(columns => [
    {
      id: 'selection',
      Header: ({ getToggleAllRowsSelectedProps, flatRows, isAllRowsSelected }) => {
        const { onChange, ...rest } = getToggleAllRowsSelectedProps();

        return (
          <Checkbox
            {...rest}
            onChange={event => {
              onChange(event);
              setRowsSelected(
                batch.id,
                flatRows.filter(row => row.depth === 1),
                !isAllRowsSelected
              );
            }}
          />
        );
      },
      Cell: ({ row }) => {
        const { onChange, ...rest } = row.getToggleRowSelectedProps();

        return (
          <Checkbox
            {...rest}
            onClick={event => event.stopPropagation()}
            onChange={event => {
              onChange(event);
              // If this is a target row, select only its subRows
              setRowsSelected(batch.id, row.depth === 0 ? row.subRows : [row], !row.isSelected);
            }}
          />
        );
      }
    },
    ...columns
  ]);
};
