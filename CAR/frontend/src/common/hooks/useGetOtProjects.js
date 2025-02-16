import { useQuery } from 'react-query';
import { getOtProjectsQueryKey } from '../api/otProjectsQueryKeys';
import { axiosGet } from '../utils/axiosFunctions';

export const useGetOtProjects = params => {
  const queryKey = getOtProjectsQueryKey(params);

  return useQuery(queryKey, async () => (await axiosGet(queryKey)).results);
};
