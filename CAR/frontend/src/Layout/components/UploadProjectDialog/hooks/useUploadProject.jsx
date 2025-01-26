import React from 'react';
import { useMutation, useQueryClient } from 'react-query';
import { axiosPost } from '../../../../common/utils/axiosFunctions';
import { addCeleryTask } from '../../../../common/stores/celeryTasksStore';
import { CloseSnackbarButton } from '../../../../common/components/CloseSnackbarButton';
import { scopes } from '../../../../common/constants/scopes';
import {
  getProjectsQueryKey,
  getProjectUploadTaskStatusQueryKey,
  uploadProjectKey
} from '../../../../common/api/projectsQueryKeys';
import { ShowProjectButton } from '../components/ShowProjectButton';
import { useTemporaryId } from '../../../../common/hooks/useTemporaryId';
import { useGlobalSnackbar } from '../../../../common/hooks/useGlobalSnackbar';

export const useUploadProject = () => {
  const queryClient = useQueryClient();

  const { generateId } = useTemporaryId();

  const { enqueueSnackbarInfo, enqueueSnackbarSuccess, enqueueSnackbarError, closeSnackbar } = useGlobalSnackbar();

  const projectsQueryKey = getProjectsQueryKey();

  return useMutation(
    ({ data }) => {
      const formData = new FormData();
      Object.entries(data).forEach(([key, value]) => {
        formData.append(key, value);
      });
      return axiosPost(uploadProjectKey(), formData);
    },
    {
      onMutate: async ({ data }) => {
        const creatingMessageId = enqueueSnackbarInfo('A project is being created...');

        // Cancel any outgoing refetches (so they don't overwrite our optimistic update)
        await queryClient.cancelQueries(projectsQueryKey);

        // Snapshot the previous value
        const previousProjects = queryClient.getQueryData(projectsQueryKey);

        // Optimistically update to the new value
        queryClient.setQueryData(projectsQueryKey, projects => {
          const newProject = {
            id: generateId(),
            init_date: null,
            name: data.project_name,
            submitterorganisation: data.submitter_organisation,
            submittername: data.submitter_name,
            proteintarget: data.protein_target,
            quotedcost: null,
            quoteurl: null
          };

          return [...projects, newProject];
        });

        // Return a context object with the snapshotted value
        return { previousProjects, creatingMessageId };
      },
      // If the mutation fails, use the context returned from onMutate to roll back
      onError: (err, vars, { previousProjects, creatingMessageId }) => {
        closeSnackbar(creatingMessageId);

        console.error(err);
        enqueueSnackbarError(err.message);

        queryClient.setQueryData(projectsQueryKey, previousProjects);

        queryClient.invalidateQueries(projectsQueryKey);
      },
      onSuccess: ({ task_id }, vars, { creatingMessageId }) => {
        addCeleryTask(task_id, {
          queryKey: getProjectUploadTaskStatusQueryKey({ task_id }),
          scope: scopes.GLOBAL,
          onSuccess: async ({ project_id }) => {
            /**
             * Load the projects again first, its used to navigate to the project in the ShowProjectButton component. The
             * reason why it's not completely optimistically updated is because some of the project fields are generated on the
             * server. Setting a temporary project as a current project could lead to bugs since the fields wouldn't match.
             * It would be possible to do full optimistic update by setting the temporary project and then replace it with
             * a new the original one once the data is loaded, but there's a time frame where the data wouldn't match.
             */
            await queryClient.refetchQueries(projectsQueryKey);

            closeSnackbar(creatingMessageId);

            enqueueSnackbarSuccess('Project has been created successfully', {
              action: key => (
                <>
                  <ShowProjectButton messageId={key} projectId={project_id} queryClient={queryClient} />
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

            queryClient.invalidateQueries(projectsQueryKey);
          }
        });
      }
    }
  );
};
