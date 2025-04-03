import React from 'react';
import { ProjectView } from '../../../Project';
import { Header } from '../Header';
import { styled } from '@mui/material/styles';
import { BatchNavigation } from '../BatchNavigation';
import { ContentBox } from '../../../common/components/ContentBox';
import { useClearStoresOnProjectChange } from './hooks/useClearStoresOnProjectChange';
import { useCurrentProjectStore } from '../../../common/stores/currentProjectStore';
import { SuspenseWithBoundary } from '../../../common/components/SuspenseWithBoundary';
import { useLayoutStore } from '../../../common/stores/layoutStore';
import { Landing } from '../Landing';
import { SmilesValidationErrorsDialog } from '../SmilesValidationErrorsDialog';

const StyledRoot = styled('div')(({ theme }) => ({
  display: 'flex',
  flexDirection: 'column',
  minHeight: '100vh'
}));

const ContentContainer = styled('div')(({ theme }) => ({
  display: 'flex',
  gap: theme.spacing(2),
  padding: theme.spacing(2)
}));

const Navigation = styled('aside')(({ theme }) => ({
  flex: '1 0 300px'
}));

const NavigationBox = styled('div')(({ theme }) => ({
  position: 'sticky',
  top: theme.spacing(2)
}));

const ProjectContainer = styled('main')(({ theme }) => ({
  flex: '1 1 100%',
  display: 'flex',
  flexDirection: 'column',
  gap: theme.spacing(2)
}));

export const Layout = () => {
  const navigationDisplayed = useLayoutStore.useNavigation();
  const projectViewDisplayed = useLayoutStore.useProjectView();
  const currentProject = useCurrentProjectStore.useCurrentProject();

  useClearStoresOnProjectChange();

  return (
    <StyledRoot>
      <Header />
      <ContentContainer>
        {currentProject ? (
          <>
            {navigationDisplayed && (
              <Navigation>
                <NavigationBox>
                  <ContentBox title={currentProject.name}>
                    <BatchNavigation />
                  </ContentBox>
                </NavigationBox>
              </Navigation>
            )}
            {projectViewDisplayed && (
              <ProjectContainer>
                <SuspenseWithBoundary>
                  <ProjectView />
                </SuspenseWithBoundary>
              </ProjectContainer>
            )}
          </>
        ) : (
          <Landing />
        )}
      </ContentContainer>
      <SmilesValidationErrorsDialog />
    </StyledRoot>
  );
};

Layout.displayName = 'Layout';
