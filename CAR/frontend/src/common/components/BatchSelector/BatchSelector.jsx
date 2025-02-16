import React from 'react';
import { ChevronRight, ExpandMore } from '@material-ui/icons';
import { makeStyles } from '@material-ui/core';
import { TreeView } from '@material-ui/lab';
import { BatchSelectorItem } from './components/BatchSelectorItem';
import { useBatchTree } from '../../hooks/useBatchTree';

const useStyles = makeStyles(theme => ({
  icon: {
    color: theme.palette.action.active
  }
}));

export const BatchSelector = ({ selectedBatchesMap, onBatchSelected }) => {
  const classes = useStyles();

  const batchTree = useBatchTree();

  const selectedBatchesIds = Object.entries(selectedBatchesMap)
    .filter(([_, value]) => value)
    .map(([key]) => String(key));

  const renderTree = node => {
    const { batch } = node;

    return (
      <BatchSelectorItem
        key={batch.id}
        batch={batch}
        selected={!!selectedBatchesMap[batch.id]}
        onSelect={onBatchSelected}
      >
        {Array.isArray(node.children) ? node.children.map(node => renderTree(node)) : null}
      </BatchSelectorItem>
    );
  };

  return (
    <TreeView
      defaultCollapseIcon={<ExpandMore className={classes.icon} />}
      defaultExpandIcon={<ChevronRight className={classes.icon} />}
      selected={selectedBatchesIds}
      disableSelection
      multiSelect
    >
      {batchTree.map(item => renderTree(item))}
    </TreeView>
  );
};
