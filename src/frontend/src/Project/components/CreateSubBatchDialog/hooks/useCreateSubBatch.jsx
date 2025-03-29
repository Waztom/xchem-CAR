import React from 'react';
import { useQueryClient, useMutation } from 'react-query';
import { createSubBatchKey, getBatchesQueryKey } from '../../../../common/api/batchesQueryKeys';
import { useTemporaryId } from '../../../../common/hooks/useTemporaryId';
import { useCurrentProjectStore } from '../../../../common/stores/currentProjectStore';
import { axiosPost } from '../../../../common/utils/axiosFunctions';
import { ShowSubBatchButton } from '../components/ShowSubBatchButton';
import { CloseSnackbarButton } from '../../../../common/components/CloseSnackbarButton';
import { useProjectSnackbar } from '../../../../common/hooks/useProjectSnackbar';
import { useBatchContext } from '../../../hooks/useBatchContext';
import { useGlobalSnackbar } from '../../../../common/hooks/useGlobalSnackbar';

export const useCreateSubBatch = () => {
  const queryClient = useQueryClient();

  const currentProject = useCurrentProjectStore.useCurrentProject();
  const batch = useBatchContext();

  const batchesQueryKey = getBatchesQueryKey({ project_id: currentProject.id });

  const { generateId } = useTemporaryId();

  const { enqueueSnackbarInfo, enqueueSnackbarSuccess, closeSnackbar } = useProjectSnackbar();
  const { enqueueSnackbarError } = useGlobalSnackbar();

  return useMutation(
    ({ batchtag, methodids }) =>
      axiosPost(createSubBatchKey(), {
        batchtag,
        methodids
      }),
    {
      onMutate: async ({ batchtag }) => {
        const creatingMessageId = enqueueSnackbarInfo('A new subbatch is being created...');

        // Cancel any outgoing refetches (so they don't overwrite our optimistic update)
        await queryClient.cancelQueries(batchesQueryKey);

        // Snapshot the previous value
        const previousBatches = queryClient.getQueryData(batchesQueryKey);

        // Temporary batchId
        const temporaryId = generateId();

        // Optimistically update to the new value
        queryClient.setQueryData(batchesQueryKey, batches => {
          const newBatch = { id: temporaryId, batch_id: batch.id, project_id: currentProject.id, batchtag: batchtag };

          return [...batches, newBatch];
        });

        // Return a context object with the snapshotted value
        return { previousBatches, temporaryId, creatingMessageId };
      },
      // If the mutation fails, use the context returned from onMutate to roll back
      onError: (err, vars, { previousBatches, creatingMessageId }) => {
        closeSnackbar(creatingMessageId);

        console.error(err);
        enqueueSnackbarError(err.message);

        queryClient.setQueryData(batchesQueryKey, previousBatches);
      },
      onSuccess: (response, vars, { temporaryId, creatingMessageId }) => {
        // Real batch ID
        const batchId = response.id;

        // Optimistically update to the new value
        queryClient.setQueryData(batchesQueryKey, batches => {
          return batches.map(batch => {
            // Temporary ID needs to be replaced with the real one
            if (batch.id === temporaryId) {
              return { ...batch, id: batchId };
            }
            return batch;
          });
        });

        closeSnackbar(creatingMessageId);

        enqueueSnackbarSuccess('The subbatch was created successfully', {
          action: key => (
            <>
              <ShowSubBatchButton queryClient={queryClient} messageId={key} batchId={batchId} />
              <CloseSnackbarButton messageId={key} />
            </>
          )
        });
      },
      // Always refetch after error or success:
      onSettled: () => {
        queryClient.invalidateQueries(batchesQueryKey);
      }
    }
  );
};
