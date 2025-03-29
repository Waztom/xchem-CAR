import create from 'zustand';
import { createSelectorHooks } from 'auto-zustand-selectors-hook';

/**
 * Stores the state of OT protocol summary dialog. You can open the dialog either from the navigation menu or through
 * a notification.
 */
const otProjectSummaryDialogStore = create(() => ({
  dialogOpen: false,
  otProjectId: null
}));

export const useOtProjectSummaryDialogStore = createSelectorHooks(otProjectSummaryDialogStore);

export const requestOtProjectSummary = otProjectId =>
  useOtProjectSummaryDialogStore.setState({ dialogOpen: true, otProjectId });

export const setOtProjectSummaryDialogOpen = dialogOpen => useOtProjectSummaryDialogStore.setState({ dialogOpen });

export const setOtProjectForSummaryDialog = otProjectId => useOtProjectSummaryDialogStore.setState({ otProjectId });
