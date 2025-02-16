import React from 'react';
import { BatchView } from '../BatchView';
import { BatchContext } from '../../context/BatchContext';
import { useBatchNavigationStore } from '../../../common/stores/batchNavigationStore';
import { useCurrentProjectStore } from '../../../common/stores/currentProjectStore';
import { useGetBatches } from '../../../common/hooks/useGetBatches';
import { ReactionDetailsDialog } from '../ReactionDetailsDialog/ReactionDetailsDialog';

export const ProjectView = () => {
  const currentProject = useCurrentProjectStore.useCurrentProject();

  const { data: batches } = useGetBatches({ project_id: currentProject.id });

  const selected = useBatchNavigationStore.useSelected();

  return (
    <>
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
    </>
  );
};
