import React, { useState } from 'react';
import ReactDOM from 'react-dom';

import Header from '../../../Layout/Header';
import Body from '../../../Body/Body';
import { ThemeProvider } from '@material-ui/core';
import Project from '../../../Project';
import { theme } from './theme';
import { QueryClient, QueryClientProvider } from 'react-query';

const queryClient = new QueryClient({
  defaultOptions: {
    queries: {
      staleTime: Infinity,
    },
  },
});

export const App = () => {
  const [projectId, setProjectId] = useState();

  function changeProject(projectid) {
    setProjectId(projectid);
  }

  function deleteProject() {
    setProjectId();
  }

  function generateProtocol() {
    setProjectId(0);
    // Need to figure out how to setshow Body vs Prot generator -> setShow in state did not work...
  }

  return (
    <QueryClientProvider client={queryClient}>
      <ThemeProvider theme={theme}>
        <Header
          changeProject={changeProject}
          deleteProject={deleteProject}
          generateProtocol={generateProtocol}
        />
        {projectId && <Project projectId={projectId} />}
        {/* <ProtocolBody ProjectID={ProjectID} key={ProjectID + 1} /> */}
      </ThemeProvider>
    </QueryClientProvider>
  );
};
