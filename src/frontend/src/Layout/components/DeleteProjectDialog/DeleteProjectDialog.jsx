import React, { useState } from 'react';
import { Typography } from '@material-ui/core';
import { useDeleteProject } from './hooks/useDeleteProject';
import { ConfirmationDialog } from '../../../common/components/ConfirmationDialog';
import {
  setDeleteProjectDialogOpen,
  setProjectForDeleteProjectDialog,
  useDeleteProjectDialogStore
} from '../../stores/deleteProjectDialogStore';

export const DeleteProjectDialog = () => {
  const { dialogOpen, project } = useDeleteProjectDialogStore();

  const [okDisabled, setOkDisabled] = useState(false);

  const { mutate: deleteProject } = useDeleteProject();

  return (
    <ConfirmationDialog
      id="delete-project-dialog"
      open={dialogOpen}
      title="Delete project"
      content={
        <Typography>
          Are you sure you want to delete project <strong>{project?.name}</strong>?
        </Typography>
      }
      onClose={() => setDeleteProjectDialogOpen(false)}
      onOk={() => {
        setOkDisabled(true);
        setDeleteProjectDialogOpen(false);
        deleteProject({ project });
      }}
      okDisabled={okDisabled}
      TransitionProps={{
        onExited: () => {
          // Prevents project name from suddenly disappearing when the dialog is closing
          setProjectForDeleteProjectDialog(null);
          setOkDisabled(false);
        }
      }}
    />
  );
};
