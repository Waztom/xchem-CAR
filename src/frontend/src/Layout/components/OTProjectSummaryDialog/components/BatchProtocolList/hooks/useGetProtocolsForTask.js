import { useMemo } from 'react';
import { useGetBatches } from '../../../../../../common/hooks/useGetBatches';
import { useGetOtBatchProtocols } from '../../../../../../common/hooks/useGetOtBatchProtocols';
import { useCurrentProjectStore } from '../../../../../../common/stores/currentProjectStore';

export const useGetProtocolsForTask = otProjectId => {
  const currentProject = useCurrentProjectStore.useCurrentProject();

  const { data: batches } = useGetBatches({ project_id: currentProject.id });
  const { data: otBatchProtocols } = useGetOtBatchProtocols({
    project_id: currentProject.id,
    otproject_id: otProjectId
  });

  return useMemo(() => {
    if (!batches || !otBatchProtocols) {
      return [];
    }

    return (
      otBatchProtocols
        .map(otBatchProtocol => ({
          otBatchProtocol,
          batch: batches.find(batch => batch.id === otBatchProtocol.batch_id)
        }))
        /**
         * If a batch in a protocol gets deleted, the optimistic update removes it from the batches data and thus the find
         * operation returns undefined. In that case, the OT Batch Protocol doesn't exists and this filter is used to bridge
         * the gap until we receive fresh data from the server where the OT Batch Protocol entry is no longer present.
         */
        .filter(({ batch }) => !!batch)
    );
  }, [batches, otBatchProtocols]);
};
