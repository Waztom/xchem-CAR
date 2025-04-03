import React from 'react';
import { AppBar, Toolbar, Typography, Button } from '@mui/material';
import { styled } from '@mui/material/styles';
import { setNavigationDisplayed, setProjectViewDisplayed, useLayoutStore } from '../../../common/stores/layoutStore';
import { setCurrentProject, useCurrentProjectStore } from '../../../common/stores/currentProjectStore';
import { LayoutSwitch } from './components/LayoutSwitch/LayoutSwitch';
import { CreateOTProjectButton } from './components/CreateOTProjectButton';
import { OTProjectHistoryButton } from './components/OTProjectHistoryButton';
import { ExportProjectButton } from './components/ExportProjectButton';

const StyledToolbar = styled(Toolbar)(({ theme }) => ({
  gap: theme.spacing(),
  '& button, a': {
    color: `${theme.palette.common.white} !important`
  }
}));

const StyledHeading = styled(Button)(({ theme }) => ({
  textTransform: 'none'
}));

export const Header = () => {
  const navigationDisplayed = useLayoutStore.useNavigation();
  const projectViewDisplayed = useLayoutStore.useProjectView();
  const currentProject = useCurrentProjectStore.useCurrentProject();

  return (
    <AppBar position="static">
      <StyledToolbar>
        <StyledHeading onClick={() => setCurrentProject(null)} disableRipple>
          <Typography variant="h6" component="h1">
            Chemist Assisted Robotics
          </Typography>
        </StyledHeading>
        {!!currentProject && (
          <>
            <CreateOTProjectButton />
            <OTProjectHistoryButton />
            <ExportProjectButton />
            <LayoutSwitch 
              checked={navigationDisplayed} 
              onChange={setNavigationDisplayed} 
              label="Navigation" 
            />
            <LayoutSwitch 
              checked={projectViewDisplayed} 
              onChange={setProjectViewDisplayed} 
              label="Project view" 
            />
          </>
        )}
      </StyledToolbar>
    </AppBar>
  );
};
