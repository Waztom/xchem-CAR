import React, { useEffect, useCallback, useRef } from 'react';
import { TreeView } from '@mui/x-tree-view/TreeView';
import { CircularProgress } from '@mui/material';
import { ChevronRight, ExpandMore } from '@mui/icons-material';
import { styled } from '@mui/material/styles';
import { useBatchTree } from '../../../common/hooks/useBatchTree';
import { NavigationItem } from './components/NavigationItem';
import { setBatchesExpanded, useBatchNavigationStore } from '../../../common/stores/batchNavigationStore';
import { DeleteSubBatchDialog } from '../DeleteSubBatchDialog';
import { SuspenseWithBoundary } from '../../../common/components/SuspenseWithBoundary';
import { getAllBatchNodeIds } from '../../../common/utils/getAllBatchNodeIds';

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
  const hasAutoExpanded = useRef(false);

  // Auto-expand all nodes on first render only
  useEffect(() => {
    if (batchTree?.length && !hasAutoExpanded.current) {
      hasAutoExpanded.current = true;
      setBatchesExpanded(getAllBatchNodeIds(batchTree));
    }
  }, [batchTree]);

  const allNodeIds = batchTree?.length ? getAllBatchNodeIds(batchTree) : [];
  const topLevelIds = batchTree?.map(n => String(n.batch.id)) || [];

  const handleNodeToggle = useCallback((event, nodeIds) => {
    // Find which top-level node was toggled
    const wasExpanded = expanded;
    const added = nodeIds.filter(id => !wasExpanded.includes(id));
    const removed = wasExpanded.filter(id => !nodeIds.includes(id));

    // If a top-level node was expanded, expand all its descendants too
    if (added.length > 0) {
      const topLevelExpanded = added.filter(id => topLevelIds.includes(id));
      if (topLevelExpanded.length > 0) {
        // Find all descendant IDs for the expanded top-level nodes
        const descendantIds = topLevelExpanded.flatMap(tlId => {
          const node = batchTree.find(n => String(n.batch.id) === tlId);
          return node ? getAllBatchNodeIds([node]) : [];
        });
        const merged = [...new Set([...nodeIds, ...descendantIds])];
        setBatchesExpanded(merged);
        return;
      }
    }

    // If a top-level node was collapsed, collapse all its descendants too
    if (removed.length > 0) {
      const topLevelCollapsed = removed.filter(id => topLevelIds.includes(id));
      if (topLevelCollapsed.length > 0) {
        const descendantIds = topLevelCollapsed.flatMap(tlId => {
          const node = batchTree.find(n => String(n.batch.id) === tlId);
          return node ? getAllBatchNodeIds([node]) : [];
        });
        const filtered = nodeIds.filter(id => !descendantIds.includes(id));
        setBatchesExpanded(filtered);
        return;
      }
    }

    setBatchesExpanded(nodeIds);
  }, [expanded, topLevelIds, batchTree]);

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
      onNodeToggle={handleNodeToggle}
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
