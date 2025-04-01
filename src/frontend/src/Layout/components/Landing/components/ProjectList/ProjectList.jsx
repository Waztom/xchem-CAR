import {
  Fab,
  List,
  ListItem,
  ListItemSecondaryAction,
  ListItemText,
  Tooltip,
  CircularProgress,
} from '@mui/material';
import { styled } from '@mui/material/styles';
import { DeleteForever } from '@mui/icons-material';
import React from 'react';
import { SuspenseWithBoundary } from '../../../../../common/components/SuspenseWithBoundary';
import { requestDeleteProject } from '../../../../stores/deleteProjectDialogStore';
import { DeleteProjectDialog } from '../../../DeleteProjectDialog';
import { setCurrentProject } from '../../../../../common/stores/currentProjectStore';
import { useGetProjects } from '../../../../../common/hooks/useGetProjects';
import { useTemporaryId } from '../../../../../common/hooks/useTemporaryId';

const ICON_SIZE = 36;

const StyledFab = styled(Fab)(({ theme }) => ({
  color: theme.palette.error.main,
  minHeight: 'unset',
  width: ICON_SIZE,
  height: ICON_SIZE,
  boxShadow: 'none !important'
}));

const StyledProgress = styled(CircularProgress)(() => ({
  display: 'block'
}));

const ProjectListContent = () => {
  const { data: projects } = useGetProjects();
  const { isTemporaryId } = useTemporaryId();

  return (
    <List>
      {projects && projects.length > 0 ? (
        projects.map(project => (
          <ListItem
            key={project.id}
            onClick={() => setCurrentProject(project)}
            disabled={isTemporaryId(project.id)}
            button
          >
            <ListItemText primary={project.name} />
            <ListItemSecondaryAction>
              {isTemporaryId(project.id) ? (
                <StyledProgress size={ICON_SIZE} />
              ) : (
                <Tooltip title="Delete project">
                  <StyledFab
                    onClick={() => requestDeleteProject(project)}
                    edge="end"
                    aria-label="delete-project"
                  >
                    <DeleteForever />
                  </StyledFab>
                </Tooltip>
              )}
            </ListItemSecondaryAction>
          </ListItem>
        ))
      ) : (
        <ListItem>
          <ListItemText primary="No projects found" />
        </ListItem>
      )}
    </List>
  );
};

const LoadingFallback = () => (
  <List>
    <ListItem>
      <StyledProgress size={ICON_SIZE} />
      <ListItemText primary="Loading projects..." />
    </ListItem>
  </List>
);

const ErrorFallback = ({ error }) => (
  <List>
    <ListItem>
      <ListItemText
        primary="Error loading projects"
        secondary={error.message}
        error
      />
    </ListItem>
  </List>
);

export const ProjectList = () => (
  <SuspenseWithBoundary
    fallback={<LoadingFallback />}
    ErrorFallbackComponent={ErrorFallback}
  >
    <ProjectListContent />
    <DeleteProjectDialog />
  </SuspenseWithBoundary>
);
