import { Typography } from '@material-ui/core';
import React, { useState } from 'react';
import { ConfirmationDialog } from '../../../common/components/ConfirmationDialog';
import {
  setDeleteSubBatchDialogOpen,
  setSubBatchForDeleteSubBatchDialog,
  useDeleteSubBatchDialogStore
} from '../../stores/deleteSubBatchDialogStore';
import { useDeleteBatch } from './hooks/useDeleteBatch';

export const DeleteSubBatchDialog = () => {
  const { dialogOpen, batch } = useDeleteSubBatchDialogStore();

  const [okDisabled, setOkDisabled] = useState(false);

  const { mutate: deleteBatch } = useDeleteBatch();

  return (
    <ConfirmationDialog
      id="delete-subbatch-dialog"
      open={dialogOpen}
      title="Delete subbatch"
      content={
        <Typography>
          Are you sure you want to delete batch <strong>{batch?.batchtag}</strong>?
        </Typography>
      }
      onClose={() => setDeleteSubBatchDialogOpen(false)}
      onOk={() => {
        setOkDisabled(true);
        deleteBatch({ batch });
        setDeleteSubBatchDialogOpen(false);
      }}
      okDisabled={okDisabled}
      TransitionProps={{
        onExited: () => {
          // Prevents batch name from suddenly disappearing when the dialog is closing
          setSubBatchForDeleteSubBatchDialog(null);
          setOkDisabled(false);
        }
      }}
    />
  );
};
