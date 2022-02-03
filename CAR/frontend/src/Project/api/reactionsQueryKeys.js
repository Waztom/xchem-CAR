export const getReactionsQueryKey = (methodId) =>
  `/api/reactions?search=${methodId}`;

export const patchReactionKey = (reactionId) => `/api/reactions/${reactionId}/`;
