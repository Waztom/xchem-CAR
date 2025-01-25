import React from 'react';
import { Checkbox, makeStyles, Tooltip, Typography } from '@material-ui/core';
import { TreeItem } from '@material-ui/lab';

const useStyles = makeStyles(theme => ({
  label: {
    display: 'flex',
    minWidth: 0,
    alignItems: 'center'
  },
  action: {
    padding: 0
  },
  leaf: {
    cursor: 'default'
  },
  name: {
    flexGrow: 1
  }
}));

export const BatchSelectorItem = ({ batch, children, selected, onSelect }) => {
  const classes = useStyles();

  return (
    <TreeItem
      classes={{ label: classes.label, content: !children.length && classes.leaf }}
      nodeId={String(batch.id)}
      label={
        <>
          <Typography className={classes.name} noWrap>
            {batch.batchtag}
          </Typography>
          <Tooltip title={selected ? 'Deselect batch' : 'Select batch'}>
            <Checkbox
              checked={selected}
              className={classes.action}
              onClick={e => e.stopPropagation()}
              onChange={(_, checked) => onSelect(batch.id, checked)}
              inputProps={{ 'aria-label': batch.batchtag }}
            />
          </Tooltip>
        </>
      }
    >
      {children}
    </TreeItem>
  );
};
