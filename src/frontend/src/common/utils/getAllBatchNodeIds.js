/**
 * Recursively collects all node IDs from a batch tree.
 * Used to expand/collapse all nodes in TreeView components.
 */
export const getAllBatchNodeIds = (nodes) => {
  const ids = [];
  const collect = (nodeList) => {
    nodeList.forEach(node => {
      ids.push(String(node.batch.id));
      if (node.children?.length) {
        collect(node.children);
      }
    });
  };
  collect(nodes);
  return ids;
};
