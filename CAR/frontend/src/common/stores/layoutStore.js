import create from 'zustand';
import { createSelectorHooks } from 'auto-zustand-selectors-hook';

/**
 * Stores displayed elements in the layout.
 */
const layoutStore = create(() => ({
  navigation: true,
  projectView: true
}));

export const useLayoutStore = createSelectorHooks(layoutStore);

export const setNavigationDisplayed = displayed => useLayoutStore.setState({ navigation: displayed });

export const setProjectViewDisplayed = displayed => useLayoutStore.setState({ projectView: displayed });

export const clearLayoutStore = () => useLayoutStore.setState({ navigation: true, projectView: true });
