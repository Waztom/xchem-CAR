import React from 'react';
import { IconButton } from '@material-ui/core';
import { useSnackbar } from 'notistack';
import { Close } from '@material-ui/icons';

export const CloseSnackbarButton = ({ messageId }) => {
  const { closeSnackbar } = useSnackbar();

  return (
    <IconButton onClick={() => closeSnackbar(messageId)} color="inherit">
      <Close />
    </IconButton>
  );
};
