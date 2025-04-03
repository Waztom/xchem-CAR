import { useQuery } from 'react-query';
import { getTargetsQueryKey } from '../api/targetsQueryKeys';
import { axiosGet } from '../utils/axiosFunctions';

export const useGetTargets = params => {
  const queryKey = getTargetsQueryKey(params);

  return useQuery(
    queryKey,
    async () => {
      try {
        const response = await axiosGet(queryKey);
        return response;
      } catch (error) {
        throw error;
      }
    },
    {
      suspense: true,
      useErrorBoundary: true,
      staleTime: 5000,
      retry: 2,
    }
  );
};
