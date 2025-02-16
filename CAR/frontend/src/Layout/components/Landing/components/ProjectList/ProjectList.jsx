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
import React from 'react';
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

export const ProjectList = () => {
  const classes = useStyles();

  const { data: projects } = useGetProjects();

  const { isTemporaryId } = useTemporaryId();

  return (
    <>
      <List>
        {!!projects?.length ? (
          projects.map(project => {
            const isProjectBeingCreated = isTemporaryId(project.id);

            return (
              <ListItem
                key={project.id}
                onClick={() => setCurrentProject(project)}
                disabled={isProjectBeingCreated}
                button
              >
                <ListItemText primary={project.name} />
                <ListItemSecondaryAction>
                  {isProjectBeingCreated ? (
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
            );
          })
        ) : (
          <ListItem>There are no projects</ListItem>
        )}
      </List>

      <DeleteProjectDialog />
    </>
  );
};
