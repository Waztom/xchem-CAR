import { useCallback } from 'react';
import { useSnackbar } from 'notistack';
import { useCurrentProjectStore } from '../stores/currentProjectStore';
import { useSnackbarBase } from './useSnackbarBase';

export const useProjectSnackbar = () => {
  const { enqueueSnackbar, closeSnackbar } = useSnackbar();

  const currentProject = useCurrentProjectStore.useCurrentProject();
  const currentProjectId = currentProject?.id;

  const enqueue = useCallback(
    (...rest) => {
      if (currentProjectId === useCurrentProjectStore.getState().currentProject?.id) {
        return enqueueSnackbar(...rest);
      }
    },
    [currentProjectId, enqueueSnackbar]
  );

  return useSnackbarBase({
    enqueueSnackbar: enqueue,
    closeSnackbar
  });
};
