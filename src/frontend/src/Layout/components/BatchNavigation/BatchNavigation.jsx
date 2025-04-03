import React from 'react';
import { TreeView } from '@mui/x-tree-view/TreeView';
import { CircularProgress } from '@mui/material';
import { ChevronRight, ExpandMore } from '@mui/icons-material';
import { styled } from '@mui/material/styles';
import { useBatchTree } from '../../../common/hooks/useBatchTree';
import { NavigationItem } from './components/NavigationItem';
import { setBatchesExpanded, useBatchNavigationStore } from '../../../common/stores/batchNavigationStore';
import { DeleteSubBatchDialog } from '../DeleteSubBatchDialog';
import { SuspenseWithBoundary } from '../../../common/components/SuspenseWithBoundary';

const StyledTreeView = styled(TreeView)(({ theme }) => ({
  '& .MuiTreeItem-iconContainer .MuiSvgIcon-root': {
    color: theme.palette.action.active
  },
  minHeight: 200
}));

const LoadingWrapper = styled('div')(({ theme }) => ({
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'center',
  padding: theme.spacing(2),
  minHeight: 200
}));

const selectedBatchesIdsSelector = state =>
  Object.entries(state.selected)
    .filter(([_, value]) => value)
    .map(([batchId]) => String(batchId));

const BatchNavigationContent = () => {
  const batchTree = useBatchTree();
  const selected = useBatchNavigationStore(selectedBatchesIdsSelector);
  const expanded = useBatchNavigationStore.useExpanded();

  const renderTree = node => {
    if (!node || !node.batch) return null;
    
    return (
      <NavigationItem key={node.batch.id} node={node}>
        {Array.isArray(node.children) && node.children.length > 0
          ? node.children.map(childNode => renderTree(childNode))
          : null}
      </NavigationItem>
    );
  };

  if (!batchTree || batchTree.length === 0) {
    return (
      <LoadingWrapper>
        <div>No batches available</div>
      </LoadingWrapper>
    );
  }

  return (
    <StyledTreeView
      defaultCollapseIcon={<ExpandMore />}
      defaultExpandIcon={<ChevronRight />}
      selected={selected}
      expanded={expanded}
      onNodeToggle={(event, nodeIds) => setBatchesExpanded(nodeIds)}
      disableSelection
      multiSelect
    >
      {batchTree.map(item => renderTree(item))}
    </StyledTreeView>
  );
};

const LoadingFallback = () => (
  <LoadingWrapper>
    <CircularProgress />
  </LoadingWrapper>
);

export const BatchNavigation = () => {
  return (
    <>
      <SuspenseWithBoundary
        fallback={<LoadingFallback />}
      >
        <BatchNavigationContent />
      </SuspenseWithBoundary>
      <DeleteSubBatchDialog />
    </>
  );
};

BatchNavigation.displayName = 'BatchNavigation';
