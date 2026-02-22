import { useQuery } from 'react-query';
import { axiosGet } from '../utils/axiosFunctions';
import { getProjectsQueryKey } from '../api/projectsQueryKeys';

export const useGetProjects = (params) => {
  const queryKey = getProjectsQueryKey(params);

  return useQuery(
    queryKey,
    async () => {
      try {
        const response = await axiosGet(queryKey);
        // DRF pagination wraps results in {count, next, previous, results}
        return response.results !== undefined ? response.results : response;
      } catch (error) {
        console.error('Error fetching projects:', error);
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
