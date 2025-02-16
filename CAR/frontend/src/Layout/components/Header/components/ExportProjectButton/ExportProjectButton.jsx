import React, { useState } from 'react';
import { Button } from '@material-ui/core';
import { ExportProjectDialog } from '../../../ExportProjectDialog';

export const ExportProjectButton = () => {
  const [dialogOpen, setDialogOpen] = useState(false);

  return (
    <>
      <Button onClick={() => setDialogOpen(true)}>Export project</Button>
      <ExportProjectDialog open={dialogOpen} onClose={() => setDialogOpen(false)} />
    </>
  );
};
