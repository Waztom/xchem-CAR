import { createContext } from 'react';

const warn = () => {
  console.warn(
    "You tried to call Legacy Suspense without SuspenseWithBoundary component. Nothing will be suspended and the errors won't be caught."
  );
};

/**
 * !!! WARNING !!!
 * useQueries doesn't support Suspense yet. This file is a part of a mechanism which fixes this issue.
 * Don't use it anywhere else apart from the places where it has already been used and don't remove it from them.
 * DON'T CHANGE IT UNLESS YOU REALLY HAVE TO - THINGS MAY BREAK!!!
 */
export const SuspendQueriesContext = createContext({
  addQueries: (id, queries) => void warn(),
  removeQueries: id => void warn()
});
