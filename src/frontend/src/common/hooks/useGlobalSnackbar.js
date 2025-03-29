import { useSnackbar } from 'notistack';
import { useSnackbarBase } from './useSnackbarBase';

export const useGlobalSnackbar = () => {
  const snackbarBase = useSnackbar();

  return useSnackbarBase(snackbarBase);
};
