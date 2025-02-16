import create from 'zustand';
import { createSelectorHooks } from 'auto-zustand-selectors-hook';

/**
 * Stores the state of OT protocol summary dialog. You can open the dialog either from the navigation menu or through
 * a notification.
 */
const reactionDetailsDialogStore = create(() => ({
  dialogOpen: false,
  reaction: null
}));

export const useReactionDetailsDialogStore = createSelectorHooks(reactionDetailsDialogStore);

export const requestReactionDetailsDialog = reaction =>
  useReactionDetailsDialogStore.setState({ dialogOpen: true, reaction });

export const setReactionDetailsDialogOpen = dialogOpen => useReactionDetailsDialogStore.setState({ dialogOpen });

export const setReactionForReactionDetailsDialog = reaction => useReactionDetailsDialogStore.setState({ reaction });
