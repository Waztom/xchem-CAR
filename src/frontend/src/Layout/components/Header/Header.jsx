import React from 'react';
import { AppBar, Toolbar, Typography, makeStyles, Button } from '@material-ui/core';
import { setNavigationDisplayed, setProjectViewDisplayed, useLayoutStore } from '../../../common/stores/layoutStore';
import { setCurrentProject, useCurrentProjectStore } from '../../../common/stores/currentProjectStore';
import { LayoutSwitch } from './components/LayoutSwitch/LayoutSwitch';
import { CreateOTProjectButton } from './components/CreateOTProjectButton';
import { OTProjectHistoryButton } from './components/OTProjectHistoryButton';
import { ExportProjectButton } from './components/ExportProjectButton';

const useStyles = makeStyles(theme => ({
  toolbar: {
    gap: theme.spacing(),
    '& button, a': {
      color: `${theme.palette.white} !important`
    }
  },
  heading: {
    textTransform: 'none'
  }
}));

export const Header = () => {
  const classes = useStyles();

  const navigationDisplayed = useLayoutStore.useNavigation();
  const projectViewDisplayed = useLayoutStore.useProjectView();

  const currentProject = useCurrentProjectStore.useCurrentProject();

  return (
    <>
      <AppBar position="static">
        <Toolbar className={classes.toolbar}>
          <Button className={classes.heading} onClick={() => setCurrentProject(null)} disableRipple>
            <Typography variant="h6" component="h1">
              Chemist Assisted Robotics
            </Typography>
          </Button>
          {!!currentProject && (
            <>
              <CreateOTProjectButton />
              <OTProjectHistoryButton />
              <ExportProjectButton />
              <LayoutSwitch checked={navigationDisplayed} onChange={setNavigationDisplayed} label="Navigation" />
              <LayoutSwitch checked={projectViewDisplayed} onChange={setProjectViewDisplayed} label="Project view" />
            </>
          )}
        </Toolbar>
      </AppBar>
    </>
  );
};
