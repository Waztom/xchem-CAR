import { useMutation, useQueryClient } from 'react-query';
import { deleteBatchKey, getBatchesQueryKey } from '../../../../common/api/batchesQueryKeys';
import { useCurrentProjectStore } from '../../../../common/stores/currentProjectStore';
import { axiosDelete } from '../../../../common/utils/axiosFunctions';
import { getOtBatchProtocolsQueryKey } from '../../../../common/api/otBatchProtocolsQueryKeys';
import { useGlobalSnackbar } from '../../../../common/hooks/useGlobalSnackbar';

export const useDeleteBatch = () => {
  const queryClient = useQueryClient();

  const currentProject = useCurrentProjectStore.useCurrentProject();

  const batchesQueryKey = getBatchesQueryKey({ project_id: currentProject.id });

  const { enqueueSnackbarError } = useGlobalSnackbar();

  return useMutation(({ batch }) => axiosDelete(deleteBatchKey(batch.id)), {
    onMutate: async ({ batch }) => {
      // Cancel any outgoing refetches (so they don't overwrite our optimistic update)
      await queryClient.cancelQueries(batchesQueryKey);

      // Snapshot the previous value
      const previousBatches = queryClient.getQueryData(batchesQueryKey);

      // Optimistically update to the new value
      queryClient.setQueryData(batchesQueryKey, oldBatches => {
        const newBatches = [...oldBatches];

        const batchIndex = oldBatches.findIndex(b => b.id === batch.id);
        newBatches.splice(batchIndex, 1);

        return newBatches;
      });

      // Return a context object with the snapshotted value
      return { previousBatches };
    },
    // If the mutation fails, use the context returned from onMutate to roll back
    onError: (err, vars, { previousBatches }) => {
      console.error(err);
      enqueueSnackbarError(err.message);

      queryClient.setQueryData(batchesQueryKey, previousBatches);
    },
    // Always refetch after error or success:
    onSettled: () => {
      queryClient.invalidateQueries(batchesQueryKey);
      // Since information about batches is in OTBatchProtocol as well, it needs to be invalidated
      queryClient.invalidateQueries(getOtBatchProtocolsQueryKey({ project_id: currentProject.id }));
    }
  });
};
