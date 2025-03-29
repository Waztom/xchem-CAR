import React, { useState } from 'react';
import { IconButton, makeStyles, Tooltip } from '@material-ui/core';
import { AddCircle } from '@material-ui/icons';
import { ContentBox } from '../../../common/components/ContentBox';
import { ProjectList } from './components/ProjectList';
import { UploadProjectDialog } from '../UploadProjectDialog';

const useStyles = makeStyles(theme => ({
  root: {
    flexGrow: 1
  }
}));

export const Landing = () => {
  const classes = useStyles();

  const [uploadProjectDialogOpen, setUploadProjectDialogOpen] = useState(false);

  return (
    <>
      <ContentBox
        title="Projects"
        PaperProps={{
          className: classes.root
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

      <UploadProjectDialog open={uploadProjectDialogOpen} onClose={() => setUploadProjectDialogOpen(false)} />
    </>
  );
};
