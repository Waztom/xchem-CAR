import { useContext } from 'react';
import { ProjectContext } from '../context/ProjectContext';

export const useProjectId = () => {
  return useContext(ProjectContext);
};
