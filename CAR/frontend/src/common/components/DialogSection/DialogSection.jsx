import { makeStyles } from '@material-ui/core';
import React from 'react';

const useStyles = makeStyles(theme => ({
  section: {
    display: 'grid',
    gap: theme.spacing()
  }
}));

export const DialogSection = ({ children }) => {
  const classes = useStyles();

  return <section className={classes.section}>{children}</section>;
};
