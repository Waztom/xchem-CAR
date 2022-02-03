import { useMemo } from 'react';
import { useQueries } from 'react-query';
import { axiosGet } from '../../../../common/utils/axiosFunctions';
import { getMethodsQueryKey } from '../../../api/methodsQueryKeys';

export const useGetMethodsForTargets = (targets) => {
  const results = useQueries(
    targets.map((target) => {
      const queryKey = getMethodsQueryKey(target.id);

      return {
        queryKey,
        queryFn: async () => {
          const methods = await axiosGet(queryKey);

          return methods.filter((method) => method.target_id === target.id);
        },
        onError: (err) => console.error(err),
      };
    })
  );

  const areAllNotFetched = results.find((result) => !result.isFetched);

  const isLoading = results.find((result) => result.isLoading);

  const methodsWithTarget = useMemo(() => {
    if (areAllNotFetched) {
      return [];
    }

    return results
      .map((result, index) => ({
        result,
        target: targets[index],
      }))
      .filter(({ result }) => result.isSuccess)
      .map(({ result, target }) =>
        result.data.map((method) => ({
          target,
          method,
        }))
      )
      .flat();
  }, [targets, results, areAllNotFetched]);

  return { isLoading, methodsWithTarget };
};
