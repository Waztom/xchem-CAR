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
        console.log('response_step3', response);
        return response;
      } catch (error) {
        console.error('Error fetching targets:', error);
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
