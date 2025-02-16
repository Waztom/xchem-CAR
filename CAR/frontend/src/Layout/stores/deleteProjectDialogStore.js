import create from 'zustand';
import { createSelectorHooks } from 'auto-zustand-selectors-hook';

/**
 * Stores the state of delete project dialog. Provides a streamlined API.
 */
const deleteProjectDialogStore = create(() => ({
  dialogOpen: false,
  project: null
}));

export const useDeleteProjectDialogStore = createSelectorHooks(deleteProjectDialogStore);

export const requestDeleteProject = project => useDeleteProjectDialogStore.setState({ dialogOpen: true, project });

export const setDeleteProjectDialogOpen = dialogOpen => useDeleteProjectDialogStore.setState({ dialogOpen });

export const setProjectForDeleteProjectDialog = project => useDeleteProjectDialogStore.setState({ project });
