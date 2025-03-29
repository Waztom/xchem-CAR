import React from 'react';
import { useMutation, useQueryClient } from 'react-query';
import { axiosPost } from '../../../../common/utils/axiosFunctions';
import { addCeleryTask } from '../../../../common/stores/celeryTasksStore';
import { CloseSnackbarButton } from '../../../../common/components/CloseSnackbarButton/CloseSnackbarButton';
import { scopes } from '../../../../common/constants/scopes';
import { useProjectSnackbar } from '../../../../common/hooks/useProjectSnackbar';
import { ShowOTProjectSummaryButton } from '../components/ShowOTProjectSummaryButton';
import { useCurrentProjectStore } from '../../../../common/stores/currentProjectStore';
import {
  createOtProjectKey,
  getOtProjectsQueryKey,
  getOtProjectTaskStatusQueryKey
} from '../../../../common/api/otProjectsQueryKeys';
import { useGlobalSnackbar } from '../../../../common/hooks/useGlobalSnackbar';

export const useCreateOTProject = () => {
  const queryClient = useQueryClient();

  const currentProject = useCurrentProjectStore.useCurrentProject();

  const { enqueueSnackbarInfo, enqueueSnackbarSuccess, closeSnackbar } = useProjectSnackbar();
  const { enqueueSnackbarError } = useGlobalSnackbar();

  const otProjectsQueryKey = getOtProjectsQueryKey({ project_id: currentProject.id });

  return useMutation(data => axiosPost(createOtProjectKey(), data), {
    onMutate: async () => {
      const creatingMessageId = enqueueSnackbarInfo('An OT protocol is being generated...');
      return { creatingMessageId };
    },
    onError: (err, vars, { creatingMessageId }) => {
      closeSnackbar(creatingMessageId);

      console.error(err);
      enqueueSnackbarError(err.message);
    },
    onSuccess: ({ task_id }, vars, { creatingMessageId }) => {
      addCeleryTask(task_id, {
        queryKey: getOtProjectTaskStatusQueryKey({ task_id }),
        scope: scopes.PROJECT,
        onSuccess: ({ otproject_id }) => {
          closeSnackbar(creatingMessageId);

          queryClient.invalidateQueries(otProjectsQueryKey);

          enqueueSnackbarSuccess('OT protocol has been generated successfully', {
            action: key => (
              <>
                <ShowOTProjectSummaryButton messageId={key} otProjectId={otproject_id} />
                <CloseSnackbarButton messageId={key} />
              </>
            )
          });
        },
        onError: err => {
          closeSnackbar(creatingMessageId);

          const message = err.traceback ?? err.message;
          console.error(message);
          enqueueSnackbarError(message);

          queryClient.invalidateQueries(otProjectsQueryKey);
        }
      });
    }
  });
};
