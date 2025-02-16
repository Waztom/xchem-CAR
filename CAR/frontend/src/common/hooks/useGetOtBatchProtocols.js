import { useQuery } from 'react-query';
import { getOtBatchProtocolsQueryKey } from '../api/otBatchProtocolsQueryKeys';
import { axiosGet } from '../utils/axiosFunctions';

export const useGetOtBatchProtocols = params => {
  const queryKey = getOtBatchProtocolsQueryKey(params);

  return useQuery(queryKey, async () => (await axiosGet(queryKey)).results);
};
