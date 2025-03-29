import { useCallback } from 'react';
import { useBatchesTableStateStore } from '../../../../../../common/stores/batchesTableStateStore';
import { useBatchContext } from '../../../../../hooks/useBatchContext';
import { useGetTableData } from '../../../../../hooks/useGetTableData';

export const useSelectedTargets = () => {
  const batch = useBatchContext();

  const targets = useGetTableData();

  const selectedTargetsIdsWithCount = useBatchesTableStateStore(
    useCallback(
      state =>
        Object.entries(state.selected[batch.id] || {})
          .filter(([_, value]) => value)
          // Method row ids are in a form targetId.methodId
          .map(([key]) => Number(key.split('.')[0]))
          // Create a map of target ids and the count of methods belonging to that id
          .reduce((targetsIds, targetId) => {
            targetsIds[targetId] = targetsIds[targetId] ? targetsIds[targetId] + 1 : 1;
            return targetsIds;
          }, {}),
      [batch.id]
    )
  );

  const selectedTargets = targets
    .filter(target => !!selectedTargetsIdsWithCount[target.id])
    .map(target => ({ target, methodsCount: selectedTargetsIdsWithCount[target.id] }));

  return selectedTargets;
};
