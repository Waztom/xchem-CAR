import { useMutation, useQueryClient } from 'react-query';
import { axiosPatch } from '../../../../common/utils/axiosFunctions';
import {
  getReactionsQueryKey,
  patchReactionKey,
} from '../../../api/reactionsQueryKeys';

export const useAdjustReactionSuccessRate = () => {
  const queryClient = useQueryClient();

  return useMutation(
    ({ reaction, successrate }) =>
      axiosPatch(patchReactionKey(reaction.id), {
        successrate,
      }),
    {
      onMutate: async ({ reaction, successrate }) => {
        const reactionQueryKey = getReactionsQueryKey(reaction.method_id);

        // Cancel any outgoing refetches (so they don't overwrite our optimistic update)
        await queryClient.cancelQueries(reactionQueryKey);

        // Snapshot the previous value
        const previousReactions = queryClient.getQueryData(reactionQueryKey);

        // Optimistically update to the new value
        queryClient.setQueryData(reactionQueryKey, (oldReactions) => {
          const newReactions = [...oldReactions];
          const newReaction = { ...reaction, successrate };

          const reactionIndex = oldReactions.findIndex((r) => r === reaction);
          newReactions.splice(reactionIndex, 1, newReaction);

          return newReactions;
        });

        // Return a context object with the snapshotted value
        return { previousReactions };
      },
      // If the mutation fails, use the context returned from onMutate to roll back
      onError: (err, { reaction }, context) => {
        console.error(err);

        queryClient.setQueryData(
          getReactionsQueryKey(reaction.method_id),
          context.previousReactions
        );
      },
      // Always refetch after error or success:
      onSettled: (data, error, { reaction }) => {
        queryClient.invalidateQueries(getReactionsQueryKey(reaction.method_id));
      },
    }
  );
};
