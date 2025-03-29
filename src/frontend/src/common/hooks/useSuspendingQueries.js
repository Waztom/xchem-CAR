import { useQueries } from 'react-query';
import { useDeepCompareMemoize } from 'use-deep-compare-effect';
// This breaks the encapsulation of files but it hides a suspense mechanism nobody should tinker with or use anywhere else
import { useSuspendQueries } from '../components/SuspenseWithBoundary/hooks/useSuspendQueries';

/**
 * useQueries doesn't support Suspense yet. Use this hook instead.
 */
export const useSuspendingQueries = queries => {
  useSuspendQueries(queries);
  return useDeepCompareMemoize(useQueries(queries));
};
