import React, { useEffect, useRef } from 'react';
import { useSnackbar } from 'notistack';
import { SnackbarButton } from '../../../../../common/components/SnackbarButton';
import {
  setBatchesExpanded,
  setBatchSelected,
  useBatchNavigationStore
} from '../../../../../common/stores/batchNavigationStore';
import { useBatchViewsRefs } from '../../../../../common/stores/batchViewsRefsStore';
import { getBatchesQueryKey } from '../../../../../common/api/batchesQueryKeys';
import { useCurrentProjectStore } from '../../../../../common/stores/currentProjectStore';

/**
 * queryClient is passed because notistack provides is defined before react-query provider, thus using useQueryClient
 * would return undefined
 */
export const ShowSubBatchButton = ({ messageId, batchId, queryClient }) => {
  const expanded = useBatchNavigationStore.useExpanded();

  // Stores the zustand subscription
  const storeSubscription = useRef();

  const currentProject = useCurrentProjectStore.useCurrentProject();

  const { closeSnackbar } = useSnackbar();

  useEffect(() => {
    return () => {
      // In case the component gets unmounted with an active subscription, unsubscribe first
      if (storeSubscription.current) {
        storeSubscription.current();
      }
    };
  }, []);

  return (
    <SnackbarButton
      onClick={() => {
        // Expand the path to the subbatch recursively
        const newExpanded = new Set(expanded);

        const batches = queryClient.getQueryData(getBatchesQueryKey({ project_id: currentProject.id }));
        let batch = batches.find(b => b.id === batchId);

        if (batch) {
          while (!!batch.batch_id) {
            // eslint-disable-next-line no-loop-func
            const parent = batches.find(b => b.id === batch.batch_id);
            newExpanded.add(String(parent.id));
            batch = parent;
          }
          setBatchesExpanded([...newExpanded]);

          // Select the batch to be displayed
          setBatchSelected(batchId, true);

          // Check if the batch has already been selected before the step above, and if so, scroll to it
          const batchRef = useBatchViewsRefs.getState().refs[batchId];
          if (!!batchRef) {
            batchRef.scrollIntoView({ behavior: 'smooth' });
          } else {
            // Subscribe to the ref changes of the specified batch and scroll to it as soon as possible
            storeSubscription.current = useBatchViewsRefs.subscribe(
              state => state.refs[batchId],
              ref => {
                if (ref) {
                  storeSubscription.current();
                  ref.scrollIntoView({ behavior: 'smooth' });
                }
              }
            );
          }
        }

        closeSnackbar(messageId);
      }}
    >
      Show subbatch
    </SnackbarButton>
  );
};
