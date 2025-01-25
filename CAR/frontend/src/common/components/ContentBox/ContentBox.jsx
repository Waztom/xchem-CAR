import React, { forwardRef } from 'react';
import { makeStyles, Paper, Typography } from '@material-ui/core';
import { SuspenseWithBoundary } from '../SuspenseWithBoundary';

const useStyles = makeStyles(theme => ({
  titleWrapper: {
    height: theme.spacing(6),
    padding: `0 ${theme.spacing(2)}px`,
    display: 'flex',
    alignItems: 'center',
    color: theme.palette.white,
    backgroundColor: theme.palette.primary.main,
    gap: theme.spacing()
  },
  title: {
    flexGrow: 1
  }
}));

export const ContentBox = forwardRef(
  ({ title, children, endAdornment, PaperProps = {}, SuspenseProps, ErrorBoundaryProps }, ref) => {
    const classes = useStyles();

    return (
      <Paper ref={ref} square {...PaperProps}>
        <div className={classes.titleWrapper}>
          <Typography className={classes.title} variant="h6" component="h2" noWrap>
            {title}
          </Typography>
          {endAdornment}
        </div>
        <SuspenseWithBoundary SuspenseProps={SuspenseProps} ErrorBoundaryProps={ErrorBoundaryProps}>
          {children}
        </SuspenseWithBoundary>
      </Paper>
    );
  }
);

ContentBox.displayName = 'ContentBox';
