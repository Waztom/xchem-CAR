import React from 'react';
import { ChevronRight, ExpandMore } from '@mui/icons-material';
import { styled } from '@mui/material/styles';
import { TreeView } from '@mui/x-tree-view/TreeView';
import { BatchSelectorItem } from './components/BatchSelectorItem';
import { useBatchTree } from '../../hooks/useBatchTree';
import { SuspenseWithBoundary } from '../../components/SuspenseWithBoundary';

const StyledTreeIcon = styled('span')(({ theme }) => ({
  color: theme.palette.action.active,
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'center'
}));

const StyledTreeView = styled(TreeView)(({ theme }) => ({
  '& .MuiTreeItem-root': {
    '& .MuiTreeItem-content': {
      padding: theme.spacing(0.5, 1)
    }
  }
}));

const BatchSelectorContent = ({ selectedBatchesMap, onBatchSelected }) => {
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
        {Array.isArray(node.children) 
          ? node.children.map(node => renderTree(node)) 
          : null}
      </BatchSelectorItem>
    );
  };

  return (
    <StyledTreeView
      defaultCollapseIcon={
        <StyledTreeIcon>
          <ExpandMore />
        </StyledTreeIcon>
      }
      defaultExpandIcon={
        <StyledTreeIcon>
          <ChevronRight />
        </StyledTreeIcon>
      }
      selected={selectedBatchesIds}
      multiSelect
      disableSelection
    >
      {batchTree.map(item => renderTree(item))}
    </StyledTreeView>
  );
};

export const BatchSelector = (props) => (
  <SuspenseWithBoundary>
    <BatchSelectorContent {...props} />
  </SuspenseWithBoundary>
);

BatchSelector.displayName = 'BatchSelector';
