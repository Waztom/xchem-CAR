import { useCallback, useMemo } from 'react';

let id = 0;

/**
 * A hook which generates temporary ids. Useful when creating new items via API and you want to display temporary items
 * in the UI while the API requests resolves and queries are invalidated.
 */
export const useTemporaryId = () => {
  const generateId = useCallback(() => `#${++id}`, []);

  const isTemporaryId = useCallback(id => id[0] === '#', []);

  return useMemo(() => ({ generateId, isTemporaryId }), [generateId, isTemporaryId]);
};
