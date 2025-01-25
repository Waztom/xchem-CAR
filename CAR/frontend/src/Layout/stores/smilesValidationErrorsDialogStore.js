import create from 'zustand';
import { createSelectorHooks } from 'auto-zustand-selectors-hook';

/**
 * Stores the state of smiles validation errors dialog. Provides a streamlined API.
 */
const smilesValidationErrorsDialogStore = create(() => ({
  dialogOpen: false,
  errors: null
}));

export const useSmilesValidationErrorsDialogStore = createSelectorHooks(smilesValidationErrorsDialogStore);

export const requestShowSmilesValidationErrors = errors =>
  useSmilesValidationErrorsDialogStore.setState({ dialogOpen: true, errors });

export const setSmilesValidationErrorsDialogOpen = dialogOpen =>
  useSmilesValidationErrorsDialogStore.setState({ dialogOpen });

export const setErrorsForSmilesValidationErrorsDialog = errors =>
  useSmilesValidationErrorsDialogStore.setState({ errors });
