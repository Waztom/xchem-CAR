import create from 'zustand';
import { createSelectorHooks } from 'auto-zustand-selectors-hook';

/**
 * Stores the state of delete subbatch dialog. Provides a streamlined API.
 */
const deleteSubBatchDialogStore = create(() => ({
  dialogOpen: false,
  batch: null
}));

export const useDeleteSubBatchDialogStore = createSelectorHooks(deleteSubBatchDialogStore);

export const requestDeleteSubBatch = batch => useDeleteSubBatchDialogStore.setState({ dialogOpen: true, batch });

export const setDeleteSubBatchDialogOpen = dialogOpen => useDeleteSubBatchDialogStore.setState({ dialogOpen });

export const setSubBatchForDeleteSubBatchDialog = batch => useDeleteSubBatchDialogStore.setState({ batch });
