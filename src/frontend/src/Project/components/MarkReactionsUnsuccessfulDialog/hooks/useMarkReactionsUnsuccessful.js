import { useMutation, useQueryClient } from 'react-query';
import { updateReactionSuccess } from '../../../../common/api/batchesQueryKeys';
import { getTargetsQueryKey } from '../../../../common/api/targetsQueryKeys';
import { useGlobalSnackbar } from '../../../../common/hooks/useGlobalSnackbar';
import { axiosPost } from '../../../../common/utils/axiosFunctions';
import { useBatchContext } from '../../../hooks/useBatchContext';

export const useMarkReactionsUnsuccessful = () => {
  const queryClient = useQueryClient();

  const batch = useBatchContext();

  const targetsQueryKey = getTargetsQueryKey({ batch_id: batch.id, fetchall: 'yes' });

  const { enqueueSnackbarError } = useGlobalSnackbar();

  return useMutation(
    ({ reaction_ids }) => {
      const formData = new FormData();
      reaction_ids.forEach(id => formData.append('reaction_ids', id));

      return axiosPost(updateReactionSuccess(), formData);
    },
    {
      onMutate: async ({ reaction_ids }) => {
        // Cancel any outgoing refetches (so they don't overwrite our optimistic update)
        await queryClient.cancelQueries(targetsQueryKey);

        // Snapshot the previous value
        const previousTargets = queryClient.getQueryData(targetsQueryKey);

        // Optimistically update to the new value
        queryClient.setQueryData(targetsQueryKey, targets => {
          return targets.map(target => ({
            ...target,
            methods: target.methods.map(method => ({
              ...method,
              reactions: method.reactions.map(reaction => {
                if (reaction_ids.includes(reaction.id)) {
                  return { ...reaction, success: false };
                }

                return reaction;
              })
            }))
          }));
        });

        // Return a context object with the snapshotted value
        return { previousTargets };
      },
      // If the mutation fails, use the context returned from onMutate to roll back
      onError: (err, _, { previousTargets }) => {
        console.error(err);
        enqueueSnackbarError(err.message);

        queryClient.setQueryData(targetsQueryKey, previousTargets);
      },
      // Always refetch after error or success:
      onSettled: () => {
        queryClient.invalidateQueries(targetsQueryKey);
      }
    }
  );
};
