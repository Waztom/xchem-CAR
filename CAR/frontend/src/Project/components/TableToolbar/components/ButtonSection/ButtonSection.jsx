import { Button, Tooltip } from '@material-ui/core';
import React, { useState } from 'react';
import { setRowsExpanded } from '../../../../../common/stores/batchesTableStateStore';
import { useBatchContext } from '../../../../hooks/useBatchContext';
import { ToolbarSection } from '../ToolbarSection';
import { CreateSubBatchDialog } from '../../../CreateSubBatchDialog';
import { MarkReactionsUnsuccessfulDialog } from '../../../MarkReactionsUnsuccessfulDialog/MarkReactionsUnsuccessfulDialog';

export const ButtonSection = ({ tableInstance, selectedMethodsIds }) => {
  const { flatRows, toggleAllRowsExpanded, preFilteredFlatRows, setAllFilters } = tableInstance;

  const batch = useBatchContext();

  const filtersApplied = flatRows.length !== preFilteredFlatRows.length;

  const createSubBatchEnabled = !!selectedMethodsIds.length;

  const [createSubBatchDialogOpen, setCreateSubBatchDialogOpen] = useState(false);
  const [markReactionsUnsuccessfulDialogOpen, setMarkReactionsUnsuccessfulDialogOpen] = useState(false);

  return (
    <>
      <ToolbarSection title="Actions">
        <Button
          fullWidth
          variant="contained"
          color="secondary"
          onClick={() => {
            toggleAllRowsExpanded(true);
            setRowsExpanded(
              batch.id,
              flatRows.filter(row => row.depth === 0),
              true
            );
          }}
        >
          Expand rows
        </Button>

        <Button
          fullWidth
          variant="contained"
          color="secondary"
          onClick={() => {
            toggleAllRowsExpanded(false);
            setRowsExpanded(
              batch.id,
              flatRows.filter(row => row.depth === 0),
              false
            );
          }}
        >
          Collapse rows
        </Button>

        <Tooltip title={!createSubBatchEnabled ? 'In order to create a SubBatch some methods have to be selected' : ''}>
          <span>
            <Button
              fullWidth
              variant="contained"
              color="secondary"
              onClick={() => {
                setCreateSubBatchDialogOpen(true);
              }}
              disabled={!createSubBatchEnabled}
            >
              Create subbatch
            </Button>
          </span>
        </Tooltip>

        <Tooltip title={!filtersApplied ? 'No filters are active' : ''}>
          <span>
            <Button
              fullWidth
              variant="contained"
              color="secondary"
              onClick={() => {
                setAllFilters([]);
              }}
              disabled={!filtersApplied}
            >
              Clear filters
            </Button>
          </span>
        </Tooltip>

        <Button
          fullWidth
          variant="contained"
          color="secondary"
          onClick={() => {
            setMarkReactionsUnsuccessfulDialogOpen(true);
          }}
        >
          Mark reactions unsuccessful
        </Button>
      </ToolbarSection>

      <CreateSubBatchDialog
        open={createSubBatchDialogOpen}
        onClose={() => setCreateSubBatchDialogOpen(false)}
        selectedMethodsIds={selectedMethodsIds}
      />
      <MarkReactionsUnsuccessfulDialog
        open={markReactionsUnsuccessfulDialogOpen}
        onClose={() => setMarkReactionsUnsuccessfulDialogOpen(false)}
      />
    </>
  );
};
