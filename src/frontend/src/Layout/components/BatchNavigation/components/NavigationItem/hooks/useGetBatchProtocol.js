import { useQuery } from 'react-query';
import { axiosGet } from '../../../../../../common/utils/axiosFunctions';

/**
 * Fetches OTBatchProtocols for a given batch and returns the most recent one.
 * Uses non-suspense so the navigation tree doesn't flash a loading state.
 */
export const useGetBatchProtocol = batchId => {
  const queryKey = ['/otbatchprotocols/', { batch_id: batchId }];

  const { data } = useQuery(queryKey, () => axiosGet(queryKey), {
    enabled: !!batchId,
    suspense: false,
    staleTime: 30000,
  });

  if (!data || data.length === 0) return null;
  return data[data.length - 1];
};
