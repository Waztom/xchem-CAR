import { List, ListItem, makeStyles, Typography } from '@material-ui/core';
import React, { memo } from 'react';
import { useSelectedTargets } from './hooks/useSelectedTargets';

const useStyles = makeStyles(theme => ({
  listItem: {
    paddingTop: 0,
    paddingBottom: 0,
    gap: theme.spacing()
  }
}));

/**
 * Renders a list of selected targets. Wrapped in a memo since it doesn't accept props and it's pointless to rerender
 * it each time a field in the Create subbatch dialog changes. Saves around 20ms on change in dev environment.
 */
export const CreateSubBatchSelectedTargetsList = memo(() => {
  const classes = useStyles();

  const selectedTargets = useSelectedTargets();

  return (
    <List disablePadding>
      {selectedTargets.map(({ target, methodsCount }) => {
        return (
          <ListItem key={target.id} className={classes.listItem}>
            <img src={target.image} width={120} height={60} alt={target.name} />
            <Typography variant="caption">
              <b>{target.name}</b>&nbsp;({methodsCount})
            </Typography>
          </ListItem>
        );
      })}
    </List>
  );
});

CreateSubBatchSelectedTargetsList.displayName = 'CreateSubBatchSelectedTargetsList';
