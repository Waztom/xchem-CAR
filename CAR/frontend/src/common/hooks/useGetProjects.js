import { useQuery } from 'react-query';
import { axiosGet } from '../utils/axiosFunctions';
import { getProjectsQueryKey } from '../api/projectsQueryKeys';

export const useGetProjects = params => {
  const queryKey = getProjectsQueryKey(params);

  return useQuery(queryKey, async () => (await axiosGet(queryKey)).results);
};
