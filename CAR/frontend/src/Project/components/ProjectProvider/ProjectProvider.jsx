import { ProjectContext } from '../../context/ProjectContext';
import { ProjectView } from '../ProjectView/ProjectView';

export const ProjectProvider = ({ projectId }) => {
  if (!projectId) {
    return null;
  }

  return (
    <ProjectContext.Provider value={projectId}>
      <ProjectView />
    </ProjectContext.Provider>
  );
};
