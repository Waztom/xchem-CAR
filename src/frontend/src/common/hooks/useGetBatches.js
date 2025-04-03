import { useQuery } from 'react-query';
import { getBatchesQueryKey } from '../api/batchesQueryKeys';
import { axiosGet } from '../utils/axiosFunctions';

export const useGetBatches = params => {
  const queryKey = getBatchesQueryKey(params);

  return useQuery(
    queryKey,
    async () => {
      try {
        const response = await axiosGet(queryKey);
        return response;
      } catch (error) {
        console.error('Error fetching batches:', error);
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
