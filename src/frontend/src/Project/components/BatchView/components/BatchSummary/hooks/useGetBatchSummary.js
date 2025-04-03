import { useGetTargets } from '../../../../../../common/hooks/useGetTargets';
import { useBatchContext } from '../../../../../hooks/useBatchContext';

export const useGetBatchSummary = () => {
  const batch = useBatchContext();

  const { data: targets } = useGetTargets({ batch_id: batch.id, fetchall: 'yes' });

  const methods = targets.map(({ methods }) => methods).flat() || [];

  return {
    targets: targets.length,
    methods: methods.length
  };
};
