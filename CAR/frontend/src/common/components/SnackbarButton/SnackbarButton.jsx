import { Button, makeStyles } from '@material-ui/core';
import React from 'react';

const useStyles = makeStyles(theme => ({
  button: {
    color: theme.palette.primary.contrastText
  }
}));

export const SnackbarButton = props => {
  const classes = useStyles();

  return <Button className={classes.button} variant="outlined" color="inherit" {...props} />;
};
