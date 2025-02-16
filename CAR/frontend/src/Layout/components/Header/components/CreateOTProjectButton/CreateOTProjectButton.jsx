import React, { useState } from 'react';
import { Button } from '@material-ui/core';
import { CreateOTProjectDialog } from '../../../CreateOTProjectDialog';

export const CreateOTProjectButton = () => {
  const [dialogOpen, setDialogOpen] = useState(false);

  return (
    <>
      <Button onClick={() => setDialogOpen(true)}>Create OT Protocol</Button>
      <CreateOTProjectDialog open={dialogOpen} onClose={() => setDialogOpen(false)} />
    </>
  );
};
