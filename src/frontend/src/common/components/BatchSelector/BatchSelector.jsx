import React, { useCallback, useEffect, useRef, useState } from 'react';
import { ChevronRight, ExpandMore } from '@mui/icons-material';
import { styled } from '@mui/material/styles';
import { TreeView } from '@mui/x-tree-view/TreeView';
import { BatchSelectorItem } from './components/BatchSelectorItem';
import { useBatchTree } from '../../hooks/useBatchTree';
import { SuspenseWithBoundary } from '../../components/SuspenseWithBoundary';
import { getAllBatchNodeIds } from '../../utils/getAllBatchNodeIds';

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
  const [expanded, setExpanded] = useState([]);
  const hasAutoExpanded = useRef(false);

  // Auto-expand all nodes on first render only
  useEffect(() => {
    if (batchTree?.length && !hasAutoExpanded.current) {
      hasAutoExpanded.current = true;
      setExpanded(getAllBatchNodeIds(batchTree));
    }
  }, [batchTree]);

  const allNodeIds = batchTree?.length ? getAllBatchNodeIds(batchTree) : [];
  const topLevelIds = batchTree?.map(n => String(n.batch.id)) || [];

  const handleNodeToggle = useCallback((event, nodeIds) => {
    const added = nodeIds.filter(id => !expanded.includes(id));
    const removed = expanded.filter(id => !nodeIds.includes(id));

    if (added.length > 0) {
      const topLevelExpanded = added.filter(id => topLevelIds.includes(id));
      if (topLevelExpanded.length > 0) {
        const descendantIds = topLevelExpanded.flatMap(tlId => {
          const node = batchTree.find(n => String(n.batch.id) === tlId);
          return node ? getAllBatchNodeIds([node]) : [];
        });
        setExpanded([...new Set([...nodeIds, ...descendantIds])]);
        return;
      }
    }

    if (removed.length > 0) {
      const topLevelCollapsed = removed.filter(id => topLevelIds.includes(id));
      if (topLevelCollapsed.length > 0) {
        const descendantIds = topLevelCollapsed.flatMap(tlId => {
          const node = batchTree.find(n => String(n.batch.id) === tlId);
          return node ? getAllBatchNodeIds([node]) : [];
        });
        setExpanded(nodeIds.filter(id => !descendantIds.includes(id)));
        return;
      }
    }

    setExpanded(nodeIds);
  }, [expanded, topLevelIds, batchTree]);

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
      expanded={expanded}
      onNodeToggle={handleNodeToggle}
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
