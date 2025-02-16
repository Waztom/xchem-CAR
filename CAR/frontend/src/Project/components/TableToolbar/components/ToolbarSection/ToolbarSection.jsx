import React from 'react';
import { makeStyles, Typography } from '@material-ui/core';

const useStyles = makeStyles(theme => ({
  root: {
    marginBottom: theme.spacing()
  },
  content: {
    display: 'grid',
    gap: theme.spacing()
  }
}));

export const ToolbarSection = ({ title, children }) => {
  const classes = useStyles();

  return (
    <section className={classes.root}>
      <Typography variant="h6" component="p">
        {title}
      </Typography>
      <div className={classes.content}>{children}</div>
    </section>
  );
};
