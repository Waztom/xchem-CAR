import React, { useState } from 'react';
import { IconButton, Tooltip } from '@mui/material';
import { styled } from '@mui/material/styles';
import { AddCircle } from '@mui/icons-material';
import { ContentBox } from '../../../common/components/ContentBox';
import { ProjectList } from './components/ProjectList';
import { UploadProjectDialog } from '../UploadProjectDialog';

const StyledRoot = styled('div')(({ theme }) => ({
  flexGrow: 1
}));

export const Landing = () => {
  const [uploadProjectDialogOpen, setUploadProjectDialogOpen] = useState(false);

  return (
    <>
      <ContentBox
        title="Projects"
        PaperProps={{
          component: StyledRoot
        }}
        endAdornment={
          <Tooltip title="Add new project">
            <IconButton color="inherit" size="small" onClick={() => setUploadProjectDialogOpen(true)}>
              <AddCircle fontSize="large" />
            </IconButton>
          </Tooltip>
        }
      >
        <ProjectList />
      </ContentBox>

      <UploadProjectDialog 
        open={uploadProjectDialogOpen} 
        onClose={() => setUploadProjectDialogOpen(false)} 
      />
    </>
  );
};
