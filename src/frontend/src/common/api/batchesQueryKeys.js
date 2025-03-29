export const getBatchesQueryKey = params => ['/batches/', params];

export const createSubBatchKey = () => '/batches/';

export const deleteBatchKey = batchId => `/batches/${batchId}/`;

export const canonicalizeSmilesKey = () => '/batches/canonicalizesmiles/';

export const getCanonicalizeSmilesTaskStatusQueryKey = params => ['/batches/gettaskstatus/', params];

export const updateReactionSuccess = () => '/batches/updatereactionsuccess/';
