import { useSnackbar } from 'notistack';
import React from 'react';
import { getProjectsQueryKey } from '../../../../../common/api/projectsQueryKeys';
import { SnackbarButton } from '../../../../../common/components/SnackbarButton';
import { setCurrentProject } from '../../../../../common/stores/currentProjectStore';

/**
 * queryClient is passed because notistack provides is defined before react-query provider, thus using useQueryClient
 * would return undefined
 */
export const ShowProjectButton = ({ messageId, projectId, queryClient }) => {
  const { closeSnackbar } = useSnackbar();

  return (
    <SnackbarButton
      onClick={() => {
        closeSnackbar(messageId);

        const projects = queryClient.getQueryData(getProjectsQueryKey());
        const project = projects?.find(project => project.id === projectId);

        // Either refetch of projects data or the project might have been deleted
        if (project) {
          setCurrentProject(project);
        }
      }}
    >
      Show project
    </SnackbarButton>
  );
};
