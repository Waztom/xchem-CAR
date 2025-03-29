import create from 'zustand';
import { createSelectorHooks } from 'auto-zustand-selectors-hook';
import { subscribeWithSelector } from 'zustand/middleware';

const batchesRefs = create(
  subscribeWithSelector(() => ({
    refs: {}
  }))
);

export const useBatchViewsRefs = createSelectorHooks(batchesRefs);

export const setBatchViewRef = (batchId, element) =>
  useBatchViewsRefs.setState(state => ({
    refs: {
      ...state.refs,
      [batchId]: element
    }
  }));

export const clearBatchViewsRefsStore = () => useBatchViewsRefs.setState({ expanded: {}, selected: {} });
