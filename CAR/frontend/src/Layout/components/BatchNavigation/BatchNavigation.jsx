import React from 'react';
import { useBatchTree } from '../../../common/hooks/useBatchTree';
import { TreeView } from '@material-ui/lab';
import { ChevronRight, ExpandMore } from '@material-ui/icons';
import { NavigationItem } from './components/NavigationItem';
import { makeStyles } from '@material-ui/core';
import { setBatchesExpanded, useBatchNavigationStore } from '../../../common/stores/batchNavigationStore';
import { DeleteSubBatchDialog } from '../DeleteSubBatchDialog';

const useStyles = makeStyles(theme => ({
  icon: {
    color: theme.palette.action.active
  }
}));

const selectedBatchesIdsSelector = state =>
  Object.entries(state.selected)
    .filter(([_, value]) => value)
    .map(([batchId]) => String(batchId));

export const BatchNavigation = () => {
  const classes = useStyles();

  const batchTree = useBatchTree();

  const selected = useBatchNavigationStore(selectedBatchesIdsSelector);
  const expanded = useBatchNavigationStore.useExpanded();

  const renderTree = node => {
    const { batch } = node;

    return (
      <NavigationItem key={batch.id} node={node}>
        {Array.isArray(node.children) ? node.children.map(node => renderTree(node)) : null}
      </NavigationItem>
    );
  };

  return (
    <>
      <TreeView
        defaultCollapseIcon={<ExpandMore className={classes.icon} />}
        defaultExpandIcon={<ChevronRight className={classes.icon} />}
        selected={selected}
        expanded={expanded}
        onNodeToggle={(event, nodeIds) => setBatchesExpanded(nodeIds)}
        disableSelection
        multiSelect
      >
        {batchTree.map(item => renderTree(item))}
      </TreeView>

      <DeleteSubBatchDialog />
    </>
  );
};
