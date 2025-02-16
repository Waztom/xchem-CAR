import { useMutation } from 'react-query';
import { useProjectSnackbar } from '../../common/hooks/useProjectSnackbar';
import { canonicalizeSmilesKey, getCanonicalizeSmilesTaskStatusQueryKey } from '../../common/api/batchesQueryKeys';
import { axiosPost } from '../../common/utils/axiosFunctions';
import { addCeleryTask } from '../../common/stores/celeryTasksStore';
import { scopes } from '../../common/constants/scopes';

export const useCanonicalizeSmiles = (onCanonicalizeStart, onCanonicalizeEnd) => {
  const { enqueueSnackbarError } = useProjectSnackbar();

  return useMutation(
    ({ data }) => {
      return axiosPost(canonicalizeSmilesKey(), data);
    },
    {
      onMutate: () => {
        onCanonicalizeStart();
      },
      onError: err => {
        console.error(err);
        enqueueSnackbarError(err.message);

        onCanonicalizeEnd();
      },
      onSuccess: ({ task_id }) => {
        addCeleryTask(task_id, {
          queryKey: getCanonicalizeSmilesTaskStatusQueryKey({ task_id }),
          scope: scopes.PROJECT,
          onSuccess: async ({ error_summary, canonicalizedsmiles }) => {
            if (!!error_summary) {
              enqueueSnackbarError(error_summary);
            }
            // If error_summary is defined, canonicalizedsmiles is undefined so no ifing is necessary
            onCanonicalizeEnd(canonicalizedsmiles);
          },
          onError: err => {
            const message = err.traceback ?? err.message;
            console.error(message);
            enqueueSnackbarError(message);

            onCanonicalizeEnd();
          }
        });
      }
    }
  );
};
