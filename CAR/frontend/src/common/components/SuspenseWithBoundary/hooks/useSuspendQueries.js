import { useContext, useLayoutEffect, useMemo } from 'react';
import { SuspendQueriesContext } from '../context/SuspendQueriesContext';

let generator = 0;

/**
 * !!! WARNING !!!
 * useQueries doesn't support Suspense yet. This file is a part of a mechanism which fixes this issue.
 * Don't use it anywhere else apart from the places where it has already been used and don't remove it from them.
 * DON'T CHANGE IT UNLESS YOU REALLY HAVE TO - THINGS MAY BREAK!!!
 */
export const useSuspendQueries = queries => {
  // Each time queries change, generate a new ID. The reason for generating a new one once queries change is that
  // we can uniquely identify each query when rendering it in QueriesSuspender component.
  // eslint-disable-next-line react-hooks/exhaustive-deps
  const id = useMemo(() => ++generator, [queries]);

  const { addQueries, removeQueries } = useContext(SuspendQueriesContext);

  // When id and queries change (they change at the same time), remove previous queries with previous ID and add
  // new queries with a new ID. useLayoutEffect is used to immediately mount the queries into the tree so the
  // the component tree immediately suspends if the queries haven't been resolved yet.
  useLayoutEffect(() => {
    addQueries(id, queries);

    return () => removeQueries(id);
  }, [id, addQueries, removeQueries, queries]);
};
