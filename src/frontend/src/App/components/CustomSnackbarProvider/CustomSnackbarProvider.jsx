import React from 'react';
import { SnackbarProvider } from 'notistack';
import { styled } from '@mui/material/styles';

const StyledSnackbarProvider = styled(SnackbarProvider)(({ theme }) => ({
  '&.SnackbarContent-root': {
    pointerEvents: 'all',
    margin: theme.spacing(6/8, 0)
  }
}));

export const CustomSnackbarProvider = ({ children }) => {
  return <StyledSnackbarProvider>{children}</StyledSnackbarProvider>;
};
