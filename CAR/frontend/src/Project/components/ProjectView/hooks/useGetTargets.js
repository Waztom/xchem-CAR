import { useQuery } from 'react-query';
import { getTargetsQueryKey } from '../../../api/targetsQueryKeys';
import { axiosGet } from '../../../../common/utils/axiosFunctions';
import { useMemo } from 'react';

export const useGetTargets = (projectId) => {
  const queryKey = getTargetsQueryKey(projectId);

  const { data, ...rest } = useQuery(queryKey, () => axiosGet(queryKey), {
    onError: (err) => console.error(err),
  });

  const targets = useMemo(() => {
    if (!data) {
      return [];
    }
    return data;
  }, [data]);

  return {
    targets,
    ...rest,
  };
};
