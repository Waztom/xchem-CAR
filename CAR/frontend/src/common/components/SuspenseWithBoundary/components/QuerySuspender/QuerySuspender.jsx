import React, { useCallback, useMemo, useState } from 'react';
import { SuspendQueriesContext } from '../../context/SuspendQueriesContext';
import { QueryWrapper } from '../QueryWrapper/QueryWrapper';

/**
 * !!! WARNING !!!
 * useQueries doesn't support Suspense yet. This file is a part of a mechanism which fixes this issue.
 * Don't use it anywhere else apart from the places where it has already been used and don't remove it from them.
 * DON'T CHANGE IT UNLESS YOU REALLY HAVE TO - THINGS MAY BREAK!!!
 */
export const QuerySuspender = ({ children }) => {
  // Store queries in a map where keys are IDs of the querying hooks
  const [queriesMap, setQueriesMap] = useState({});

  const addQueries = useCallback((id, queries) => {
    setQueriesMap(prevQueriesMap => ({
      ...prevQueriesMap,
      [id]: queries
    }));
  }, []);

  const removeQueries = useCallback(id => {
    setQueriesMap(prevQueriesMap => {
      const newQueriesMap = { ...prevQueriesMap };
      delete newQueriesMap[id];
      return newQueriesMap;
    });
  }, []);

  const value = useMemo(
    () => ({
      addQueries,
      removeQueries
    }),
    [addQueries, removeQueries]
  );

  return (
    <SuspendQueriesContext.Provider value={value}>
      {children}
      {Object.entries(queriesMap).map(([id, queries]) =>
        queries.map((query, index) => <QueryWrapper key={`${id}-${index}`} query={query} />)
      )}
    </SuspendQueriesContext.Provider>
  );
};
