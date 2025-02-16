import React from 'react';
import { SnackbarProvider } from 'notistack';
import { makeStyles } from '@material-ui/core';

const useStyles = makeStyles(theme => ({
  root: {
    pointerEvents: 'all',
    margin: theme.spacing(6 / 8, 0)
  }
}));

export const CustomSnackbarProvider = ({ children }) => {
  const classes = useStyles();

  return <SnackbarProvider className={classes.root}>{children}</SnackbarProvider>;
};
