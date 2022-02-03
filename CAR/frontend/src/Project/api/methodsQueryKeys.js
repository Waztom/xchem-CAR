export const getMethodsQueryKey = (targetId) =>
  `/api/methods?search=${targetId}`;

export const patchMethodsKey = (methodId) => `/api/methods/${methodId}/`;
