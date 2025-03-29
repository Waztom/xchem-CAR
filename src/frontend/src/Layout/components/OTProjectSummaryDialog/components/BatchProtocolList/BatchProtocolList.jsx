import { Button, List, ListItem, ListItemSecondaryAction, ListItemText } from '@material-ui/core';
import React from 'react';
import { useGetProtocolsForTask } from './hooks/useGetProtocolsForTask';

export const BatchProtocolList = ({ otProjectId }) => {
  const batchesWithOtProjects = useGetProtocolsForTask(otProjectId);

  return (
    <List>
      {!!batchesWithOtProjects.length ? (
        batchesWithOtProjects.map(({ batch, otBatchProtocol }) => {
          return (
            <ListItem key={batch.id}>
              <ListItemText primary={batch.batchtag} />
              <ListItemSecondaryAction>
                <Button color="primary" variant="contained" href={otBatchProtocol.zipfile} download>
                  Download protocol
                </Button>
              </ListItemSecondaryAction>
            </ListItem>
          );
        })
      ) : (
        <ListItem>
          <ListItemText primary="None of the currently existing batches are left in the OT Protocol" />
        </ListItem>
      )}
    </List>
  );
};
