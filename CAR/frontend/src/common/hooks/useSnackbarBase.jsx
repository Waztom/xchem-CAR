import React, { useMemo } from 'react';
import { useCallback } from 'react';
import { CloseSnackbarButton } from '../components/CloseSnackbarButton/CloseSnackbarButton';

/**
 * Serves as a base for useProjectSnackbar and useGlobalSnackbar. This file should not be used directly. Use
 * useProjectSnackbar or useGlobalSnackbar instead.
 */
export const useSnackbarBase = snackbarBase => {
  const { enqueueSnackbar, closeSnackbar } = snackbarBase;

  const enqueueSnackbarInfo = useCallback(
    (message, options = {}) => {
      return enqueueSnackbar(message, {
        variant: 'info',
        ...options
      });
    },
    [enqueueSnackbar]
  );

  const enqueueSnackbarSuccess = useCallback(
    (message, options = {}) => {
      return enqueueSnackbar(message, {
        variant: 'success',
        autoHideDuration: null,
        action: key => <CloseSnackbarButton messageId={key} />,
        ...options
      });
    },
    [enqueueSnackbar]
  );

  const enqueueSnackbarError = useCallback(
    (message, options = {}) => {
      return enqueueSnackbar(message, {
        variant: 'error',
        autoHideDuration: null,
        action: key => <CloseSnackbarButton messageId={key} />,
        ...options
      });
    },
    [enqueueSnackbar]
  );

  return useMemo(
    () => ({
      enqueueSnackbar,
      enqueueSnackbarInfo,
      enqueueSnackbarSuccess,
      enqueueSnackbarError,
      closeSnackbar
    }),
    [enqueueSnackbar, enqueueSnackbarInfo, enqueueSnackbarSuccess, enqueueSnackbarError, closeSnackbar]
  );
};
