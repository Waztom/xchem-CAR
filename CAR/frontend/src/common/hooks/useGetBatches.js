import { useQuery } from 'react-query';
import { getBatchesQueryKey } from '../api/batchesQueryKeys';
import { axiosGet } from '../utils/axiosFunctions';

export const useGetBatches = params => {
  const queryKey = getBatchesQueryKey(params);

  return useQuery(queryKey, async () => (await axiosGet(queryKey)).results);
};
