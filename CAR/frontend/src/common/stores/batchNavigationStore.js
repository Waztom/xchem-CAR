import create from 'zustand';
import { createSelectorHooks } from 'auto-zustand-selectors-hook';

const batchNavigationStore = create(() => ({ selected: {}, expanded: [] }));

export const useBatchNavigationStore = createSelectorHooks(batchNavigationStore);

export const setBatchSelected = (batchId, selected) =>
  useBatchNavigationStore.setState(state => ({ selected: { ...state.selected, [batchId]: selected } }));

export const setBatchesExpanded = expanded => useBatchNavigationStore.setState(state => ({ expanded }));

export const clearBatchNavigationStore = () => useBatchNavigationStore.setState({ selected: {}, expanded: [] });
