import React from 'react';
import { makeStyles, Tooltip, Typography } from '@material-ui/core';
import { useGetBatchSummary } from './hooks/useGetBatchSummary';
import { IconComponent } from '../../../../../common/components/IconComponent';
import { FaFlask } from 'react-icons/fa';
import { FindInPage } from '@material-ui/icons';

const useStyles = makeStyles(theme => ({
  root: {
    display: 'flex'
  },
  categoryInfo: {
    display: 'grid',
    gridTemplateColumns: 'repeat(2, auto)',
    alignItems: 'center',
    gap: theme.spacing(1 / 2)
  }
}));

export const BatchSummary = () => {
  const classes = useStyles();

  const { targets, methods } = useGetBatchSummary();

  return (
    <div className={classes.root}>
      <Tooltip title={`There are ${targets} targets in total`}>
        <div className={classes.categoryInfo}>
          <Typography>{targets}</Typography>
          <IconComponent Component={FindInPage} />
        </div>
      </Tooltip>

      <Tooltip title={`There are ${methods} methods in total`}>
        <div className={classes.categoryInfo}>
          <Typography>{methods}</Typography>
          <IconComponent Component={FaFlask} />
        </div>
      </Tooltip>
    </div>
  );
};
