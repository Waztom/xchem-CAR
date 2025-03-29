/**
 * Defines scopes which can be used for data, actions, or other stuff in the app. Currently used to distinguish
 * between Celery tasks scope (e.g. Global Celery tasks shouldnt be cleared when a project is changed)
 */
export const scopes = {
  GLOBAL: 'global',
  PROJECT: 'project'
};
