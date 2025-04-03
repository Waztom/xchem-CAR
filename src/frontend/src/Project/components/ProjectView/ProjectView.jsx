import React, { Suspense } from 'react';
import { BatchView } from '../BatchView';
import { BatchContext } from '../../context/BatchContext';
import { useBatchNavigationStore } from '../../../common/stores/batchNavigationStore';
import { useCurrentProjectStore } from '../../../common/stores/currentProjectStore';
import { useGetBatches } from '../../../common/hooks/useGetBatches';
import { ReactionDetailsDialog } from '../ReactionDetailsDialog/ReactionDetailsDialog';
import { List, ListItem, ListItemText, CircularProgress } from '@mui/material';

export const ProjectView = () => {
  const currentProject = useCurrentProjectStore.useCurrentProject();
  const { data: batches } = useGetBatches({ project_id: currentProject.id });
  const selected = useBatchNavigationStore.useSelected();

  return (
    <Suspense
      fallback={
        <List>
          <ListItem>
            <CircularProgress size={36} />
            <ListItemText primary="Loading batches..." />
          </ListItem>
        </List>
      }
    >
      {batches
        ?.filter(batch => selected[batch.id])
        .map(batch => {
          return (
            <BatchContext.Provider key={batch.id} value={batch}>
              <BatchView />
            </BatchContext.Provider>
          );
        })}

      <ReactionDetailsDialog />
    </Suspense>
  );
};
