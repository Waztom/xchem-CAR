import { useMemo } from 'react';
import { useQueries } from 'react-query';
import { axiosGet } from '../../../../common/utils/axiosFunctions';
import { getReactionsQueryKey } from '../../../api/reactionsQueryKeys';

export const useGetMethodsReactions = (methodsWithTarget) => {
  const results = useQueries(
    methodsWithTarget.map((methodWithTarget) => {
      const queryKey = getReactionsQueryKey(methodWithTarget.method.id);

      return {
        queryKey,
        queryFn: async () => {
          const reactions = await axiosGet(queryKey);
          return reactions.sort(
            (reactionA, reactionB) => reactionA.id - reactionB.id
          );
        },
        onError: (err) => console.error(err),
      };
    })
  );

  const areAllNotFetched = results.find((result) => !result.isFetched);

  const isLoading = results.find((result) => result.isLoading);

  const methodsData = useMemo(() => {
    if (areAllNotFetched) {
      return [];
    }

    return results
      .map((result, index) => ({
        result,
        methodWithTarget: methodsWithTarget[index],
      }))
      .filter(({ result }) => result.isSuccess)
      .map(({ result, methodWithTarget }) => ({
        ...methodWithTarget,
        reactions: result.data,
      }));
  }, [methodsWithTarget, results, areAllNotFetched]);

  return { isLoading, methodsData };
};
