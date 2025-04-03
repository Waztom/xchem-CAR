export const getProjectsQueryKey = params => ['/projects/', params];

export const uploadProjectKey = () => '/projects/createproject/';

export const deleteProjectKey = projectId => `/projects/${projectId}/`;

export const getProjectUploadTaskStatusQueryKey = params => ['/projects/gettaskstatus/', params];
