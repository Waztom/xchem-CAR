import create from 'zustand';
import { subscribeWithSelector } from 'zustand/middleware';
import { createSelectorHooks } from 'auto-zustand-selectors-hook';

/**
 * Stores currently selected project
 */
const currentProjectStore = create(
  subscribeWithSelector(() => ({
    currentProject: null
  }))
);

export const useCurrentProjectStore = createSelectorHooks(currentProjectStore);

export const setCurrentProject = project => useCurrentProjectStore.setState({ currentProject: project });
