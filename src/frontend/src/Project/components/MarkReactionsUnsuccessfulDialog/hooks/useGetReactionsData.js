import { useMemo } from 'react';
import { useGetTargets } from '../../../../common/hooks/useGetTargets';
import { useBatchContext } from '../../../hooks/useBatchContext';

/**
 * Apart from the ids which are used mainly for the autocomplete options, an object where keys are reaction IDs
 * and values are reactions is returned is used for validating the uploaded IDs from a file as well as adding
 * more info to the autocomplete items.
 */
export const useGetReactionsData = () => {
  const batch = useBatchContext();

  const { data: targets } = useGetTargets({ batch_id: batch.id, fetchall: 'yes' });

  return useMemo(() => {
    if (!targets) {
      return {
        reactionsMap: {},
        reactionsIds: []
      };
    }

    const reactions = targets.map(({ methods }) => methods.map(({ reactions }) => reactions)).flat(2);

    return {
      reactionsMap: Object.fromEntries(reactions.map(reaction => [String(reaction.id), reaction])),
      reactionsIds: reactions.map(reaction => String(reaction.id))
    };
  }, [targets]);
};
