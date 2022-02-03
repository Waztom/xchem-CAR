import { useMutation, useQueryClient } from 'react-query';
import { axiosPatch } from '../../../../common/utils/axiosFunctions';
import {
  getMethodsQueryKey,
  patchMethodsKey,
} from '../../../api/methodsQueryKeys';

export const useSynthesiseMethod = () => {
  const queryClient = useQueryClient();

  return useMutation(
    ({ method, synthesise }) =>
      axiosPatch(patchMethodsKey(method.id), { synthesise }),
    {
      onMutate: async ({ method, synthesise }) => {
        const methodQueryKey = getMethodsQueryKey(method.target_id);

        // Cancel any outgoing refetches (so they don't overwrite our optimistic update)
        await queryClient.cancelQueries(methodQueryKey);

        // Snapshot the previous value
        const previousMethods = queryClient.getQueryData(methodQueryKey);

        // Optimistically update to the new value
        queryClient.setQueryData(methodQueryKey, (oldMethods) => {
          const newMethods = [...oldMethods];
          const newMethod = { ...method, synthesise };

          const methodIndex = oldMethods.findIndex((m) => m === method);
          newMethods.splice(methodIndex, 1, newMethod);

          return newMethods;
        });

        // Return a context object with the snapshotted value
        return { previousMethods };
      },
      // If the mutation fails, use the context returned from onMutate to roll back
      onError: (err, { method }, context) => {
        console.error(err);

        queryClient.setQueryData(
          getMethodsQueryKey(method.target_id),
          context.previousMethods
        );
      },
      // Always refetch after error or success:
      onSettled: (data, error, { method }) => {
        queryClient.invalidateQueries(getMethodsQueryKey(method.target_id));
      },
    }
  );
};
