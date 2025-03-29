import React from 'react';
import { List, ListItem, makeStyles, Typography } from '@material-ui/core';
import { Alert } from '@material-ui/lab';
import { useGetIncompatibleTargets } from './hooks/useGetIncompatibleTargets';

const useStyles = makeStyles(theme => ({
  listItem: {
    paddingTop: 0,
    paddingBottom: 0,
    gap: theme.spacing()
  },
  image: {
    mixBlendMode: 'multiply'
  }
}));

export const OTWarningSection = ({ selectedBatchesMap }) => {
  const classes = useStyles();

  const incompatibleTargets = useGetIncompatibleTargets(selectedBatchesMap);

  return (
    <>
      {incompatibleTargets.map(({ batch, targets }) => (
        <Alert key={batch.id} severity="warning" variant="outlined">
          These targets from batch <strong>{batch.batchtag}</strong> contain only methods which can't be executed on
          OpenTrons:
          <List disablePadding>
            {targets.map(target => {
              return (
                <ListItem key={target.id} className={classes.listItem}>
                  <img className={classes.image} src={target.image} width={120} height={60} alt={target.name} />
                  <Typography variant="caption">
                    <b>{target.name}</b>
                  </Typography>
                </ListItem>
              );
            })}
          </List>
        </Alert>
      ))}
    </>
  );
};
