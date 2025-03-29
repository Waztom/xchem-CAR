import { useQuery } from 'react-query';
import { getOtProjectsQueryKey } from '../api/otProjectsQueryKeys';
import { axiosGet } from '../utils/axiosFunctions';

export const useGetOtProjects = params => {
  const queryKey = getOtProjectsQueryKey(params);

  return useQuery(
    queryKey,
    async () => {
      try {
        const response = await axiosGet(queryKey);
        return response;
      } catch (error) {
        console.error('Error fetching OT projects:', error);
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
