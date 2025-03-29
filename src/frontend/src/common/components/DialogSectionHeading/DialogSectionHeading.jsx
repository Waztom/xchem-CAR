import React from 'react';
import { makeStyles, Typography } from '@material-ui/core';

const useStyles = makeStyles(theme => ({
  heading: {
    fontSize: '0.9rem',
    fontWeight: 500
  }
}));

export const DialogSectionHeading = ({ children }) => {
  const classes = useStyles();

  return <Typography className={classes.heading}>{children}</Typography>;
};
