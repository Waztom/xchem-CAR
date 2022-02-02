import { useProjectId } from '../../../hooks/useProjectId';
import { useQuery } from 'react-query';
import { getTargetsQueryKey } from '../../../api/targetsQueryKeys';
import { axiosGet } from '../../../../common/utils/axiosFunctions';

export const useGetTargets = () => {
  const projectId = useProjectId();

  const queryKey = getTargetsQueryKey(projectId);

  return useQuery(queryKey, () => axiosGet(queryKey), {
    initialData: [],
    onError: (err) => console.error(err),
  });
};
