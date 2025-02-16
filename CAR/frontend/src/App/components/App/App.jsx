import { createTheme, ThemeProvider } from '@material-ui/core';
import React from 'react';
import { QueryClient, QueryClientProvider } from 'react-query';
import { getTheme } from '../../../theme';
import Layout from '../../../Layout';
import { CeleryTasksChecker } from '../CeleryTasksChecker/CeleryTasksChecker';
import { CustomSnackbarProvider } from '../CustomSnackbarProvider';

const theme = createTheme({
  ...getTheme(),
  overrides: {
    MuiAccordionSummary: {
      root: {
        '&$expanded': {
          minHeight: 48
        }
      },
      content: {
        margin: 0,
        '&$expanded': {
          margin: 0
        }
      }
    }
  }
});

const queryClient = new QueryClient({
  defaultOptions: {
    queries: {
      staleTime: Infinity,
      suspense: true
    }
  }
});

export const App = () => {
  return (
    <ThemeProvider theme={theme}>
      <CustomSnackbarProvider>
        <QueryClientProvider client={queryClient}>
          <Layout />
          <CeleryTasksChecker />
        </QueryClientProvider>
      </CustomSnackbarProvider>
    </ThemeProvider>
  );
};
