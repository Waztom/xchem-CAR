import React from 'react';
import { IconButton } from '@mui/material';
import { useSnackbar } from 'notistack';
import { Close } from '@mui/icons-material';

export const CloseSnackbarButton = ({ messageId }) => {
  const { closeSnackbar } = useSnackbar();

  return (
    <IconButton onClick={() => closeSnackbar(messageId)} color="inherit" size="large">
      <Close />
    </IconButton>
  );
};
