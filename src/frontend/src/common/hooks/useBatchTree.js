import { useGetBatches } from './useGetBatches';
import { useCurrentProjectStore } from '../stores/currentProjectStore';

export const useBatchTree = () => {
  const currentProject = useCurrentProjectStore.useCurrentProject();

  const { data: batches } = useGetBatches({ project_id: currentProject.id });

  if (!batches) {
    return [];
  }

  const mappedNodes = {};
  const parentNodes = [];
  batches
    .sort((a, b) => a.id - b.id)
    .forEach(batch => {
      const node = { batch, children: [] };
      mappedNodes[batch.id] = node;

      if (batch.batch_id) {
        mappedNodes[batch.batch_id].children.push(node);
      } else {
        parentNodes.push(node);
      }
    });

  return parentNodes;
};
