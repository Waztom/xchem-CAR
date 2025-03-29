import {
  Fab,
  List,
  ListItem,
  ListItemSecondaryAction,
  ListItemText,
  makeStyles,
  Tooltip,
  CircularProgress
} from '@material-ui/core';
import { DeleteForever } from '@material-ui/icons';
import React, { Suspense } from 'react';
import { ErrorBoundary } from 'react-error-boundary';
import { requestDeleteProject } from '../../../../stores/deleteProjectDialogStore';
import { DeleteProjectDialog } from '../../../DeleteProjectDialog';
import { setCurrentProject } from '../../../../../common/stores/currentProjectStore';
import { useGetProjects } from '../../../../../common/hooks/useGetProjects';
import { useTemporaryId } from '../../../../../common/hooks/useTemporaryId';

const ICON_SIZE = 36;

const useStyles = makeStyles(theme => ({
  deleteButton: {
    color: theme.palette.error.main,
    minHeight: 'unset',
    width: ICON_SIZE,
    height: ICON_SIZE,
    boxShadow: 'none !important'
  },
  progress: {
    display: 'block'
  }
}));

const ProjectListContent = () => {
  const { data: projects, isLoading, isError, error } = useGetProjects();
  const classes = useStyles();
  const { isTemporaryId } = useTemporaryId();

  if (isLoading) {
    return (
      <List>
        <ListItem>
          <CircularProgress className={classes.progress} size={ICON_SIZE} />
          <ListItemText primary="Loading projects..." />
        </ListItem>
      </List>
    );
  }

  if (isError) {
    throw error;
  }

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
                <CircularProgress className={classes.progress} size={ICON_SIZE} />
              ) : (
                <Tooltip title="Delete project">
                  <Fab
                    className={classes.deleteButton}
                    onClick={() => requestDeleteProject(project)}
                    edge="end"
                    aria-label="delete-project"
                  >
                    <DeleteForever />
                  </Fab>
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

const ErrorFallback = ({ error }) => (
  <List>
    <ListItem>
      <ListItemText
        primary="Error loading projects"
        secondary={error.message}
        className="error-message"
      />
    </ListItem>
  </List>
);

export const ProjectList = () => {
  return (
    <ErrorBoundary FallbackComponent={ErrorFallback}>
      <Suspense
        fallback={
          <List>
            <ListItem>
              <CircularProgress size={ICON_SIZE} />
              <ListItemText primary="Loading projects..." />
            </ListItem>
          </List>
        }
      >
        <ProjectListContent />
        <DeleteProjectDialog />
      </Suspense>
    </ErrorBoundary>
  );
};
