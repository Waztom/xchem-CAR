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

          return {
            ...methodWithTarget,
            reactions,
          };
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
      .filter((result) => result.isSuccess)
      .map((result) => result.data);
  }, [results, areAllNotFetched]);

  return { isLoading, methodsData };
};
