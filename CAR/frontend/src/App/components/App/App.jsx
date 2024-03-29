import React, { useState } from 'react';
import ReactDOM from 'react-dom';

import Header from '../../../Layout/Header';
import Body from '../../../Body/Body';
import { ThemeProvider } from '@material-ui/core';
import { ProjectView } from '../../../Project';
import { theme } from './theme';

export const App = () => {
  const [projectId, setProjectId] = useState(0);

  function changeProject(projectid) {
    setProjectId(projectid);
  }

  function deleteProject() {
    setProjectId(0);
  }

  function generateProtocol() {
    setProjectId(0);
    // Need to figure out how to setshow Body vs Prot generator -> setShow in state did not work...
  }

  return (
    <ThemeProvider theme={theme}>
      <Header
        changeProject={changeProject}
        deleteProject={deleteProject}
        generateProtocol={generateProtocol}
      />
      <ProjectView projectId={projectId} />
      {/* <ProtocolBody ProjectID={ProjectID} key={ProjectID + 1} /> */}
    </ThemeProvider>
  );
};
