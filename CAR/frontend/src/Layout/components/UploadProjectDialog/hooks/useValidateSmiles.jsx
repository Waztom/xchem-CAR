import React from 'react';
import { useMutation } from 'react-query';
import { axiosPost } from '../../../../common/utils/axiosFunctions';
import { addCeleryTask } from '../../../../common/stores/celeryTasksStore';
import { scopes } from '../../../../common/constants/scopes';
import { getProjectUploadTaskStatusQueryKey, uploadProjectKey } from '../../../../common/api/projectsQueryKeys';
import { useGlobalSnackbar } from '../../../../common/hooks/useGlobalSnackbar';
import { CloseSnackbarButton } from '../../../../common/components/CloseSnackbarButton';
import { ShowValidationErrorsButton } from '../components/ShowValidationErrorsButton/ShowValidationErrorsButton';

export const useValidateSmiles = () => {
  const { enqueueSnackbarInfo, enqueueSnackbarSuccess, enqueueSnackbarError, closeSnackbar } = useGlobalSnackbar();

  return useMutation(
    ({ data }) => {
      const formData = new FormData();
      Object.entries(data).forEach(([key, value]) => {
        formData.append(key, value);
      });
      return axiosPost(uploadProjectKey(), formData);
    },
    {
      onMutate: async () => {
        const creatingMessageId = enqueueSnackbarInfo('The smiles file is being validated...');

        return { creatingMessageId };
      },
      onError: (err, vars, { creatingMessageId }) => {
        closeSnackbar(creatingMessageId);

        console.error(err);
        enqueueSnackbarError(err.message);
      },
      onSuccess: ({ task_id }, vars, { creatingMessageId }) => {
        addCeleryTask(task_id, {
          queryKey: getProjectUploadTaskStatusQueryKey({ task_id }),
          scope: scopes.GLOBAL,
          onSuccess: async ({ validated, validation_errors }) => {
            closeSnackbar(creatingMessageId);

            if (validated) {
              enqueueSnackbarSuccess('The smiles file has been successfully validated');
            } else {
              enqueueSnackbarError('The smiles file validation failed', {
                action: key => (
                  <>
                    <ShowValidationErrorsButton messageId={key} errors={JSON.parse(validation_errors)} />
                    <CloseSnackbarButton messageId={key} />
                  </>
                )
              });
            }
          },
          onError: err => {
            closeSnackbar(creatingMessageId);

            const message = err.traceback ?? err.message;
            console.error(message);
            enqueueSnackbarError(message);
          }
        });
      }
    }
  );
};
