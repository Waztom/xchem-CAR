import { useMemo } from 'react';
import { useGetBatches } from '../../../../../../common/hooks/useGetBatches';
import { useCurrentProjectStore } from '../../../../../../common/stores/currentProjectStore';
import { useSuspendingQueries } from '../../../../../../common/hooks/useSuspendingQueries';
import { getTargetsQueryKey } from '../../../../../../common/api/targetsQueryKeys';
import { axiosGet } from '../../../../../../common/utils/axiosFunctions';
import { use } from 'react';

export const useGetIncompatibleTargets = selectedBatchesMap => {
  const currentProject = useCurrentProjectStore.useCurrentProject();
  const { data: batches } = useGetBatches({ project_id: currentProject.id });

  const allTargetsResponses = useSuspendingQueries(
    useMemo(() => {
      if (!batches?.length) return [];
      return batches.map(batch => ({
        queryKey: getTargetsQueryKey({ batch_id: batch.id, fetchall: 'yes' }),
        queryFn: async () => {
          const response = await axiosGet(getTargetsQueryKey({ batch_id: batch.id, fetchall: 'yes' }));
          return response?.results || [];
        },
        suspense: true, 
        useErrorBoundary: true, 
        staleTime: 5000,
        retry: 2 
      }));
    }, [batches])
  );

  const allSelectedTargets = allTargetsResponses
    // Pair batches with their appropriate target response
    .map((response, index) => ({ batch: batches[index], response }))
    // Get only the ones selected
    .filter((_, index) => !!selectedBatchesMap[batches[index].id])
    // Get only target responses which has finished successfully
    .filter(({ response }) => response.isSuccess)
    // Pair the batch with the data about targets
    .map(({ batch, response }) => ({ batch, targets: response.data }));

  const incompatibleTargets = allSelectedTargets
    .map(({ batch, targets }) => ({
      batch,
      // Get only targets whose every method is not compatible with OT
      targets: targets.filter(target => target.methods.every(method => !method.otchem))
    }))
    .filter(({ targets }) => !!targets.length);

  return incompatibleTargets;
};
