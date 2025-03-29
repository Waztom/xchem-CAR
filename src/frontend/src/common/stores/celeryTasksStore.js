import create from 'zustand';
import { createSelectorHooks } from 'auto-zustand-selectors-hook';
import { scopes } from '../constants/scopes';

const defaultTaskOptions = {
  queryKey: [],
  scope: scopes.GLOBAL,
  onSuccess: () => {},
  onError: () => {}
};

/**
 * Stores celery tasks which should be polled for completion.
 */
const celeryTasksStore = create(() => ({
  tasks: {}
}));

export const useCeleryTasksStore = createSelectorHooks(celeryTasksStore);

export const addCeleryTask = (taskId, options) =>
  useCeleryTasksStore.setState(state => ({
    tasks: {
      ...state.tasks,
      [taskId]: { ...defaultTaskOptions, ...options, id: taskId }
    }
  }));

export const removeCeleryTask = taskId =>
  useCeleryTasksStore.setState(state => {
    const tasks = { ...state.tasks };
    delete tasks[taskId];
    return { tasks };
  });

export const removeCeleryTasksByScope = scope =>
  useCeleryTasksStore.setState(state => ({
    tasks: Object.fromEntries(Object.entries(state.tasks).filter(([_, value]) => value.scope !== scope))
  }));

export const clearCeleryTasksStore = () => useCeleryTasksStore.setState({ tasks: {} });
