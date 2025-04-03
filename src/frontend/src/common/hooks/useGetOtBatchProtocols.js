import { useQuery } from 'react-query';
import { getOtBatchProtocolsQueryKey } from '../api/otBatchProtocolsQueryKeys';
import { axiosGet } from '../utils/axiosFunctions';

export const useGetOtBatchProtocols = params => {
  const queryKey = getOtBatchProtocolsQueryKey(params);

  return useQuery(
    queryKey,
    async () => {
      try {
        const response = await axiosGet(queryKey);
        return response;
      } catch (error) {
        console.error('Error fetching OT Batch protocols:', error);
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
