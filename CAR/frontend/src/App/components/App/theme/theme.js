import { createTheme } from '@material-ui/core';
import palette from './palette';

export const theme = createTheme({
  palette,
  zIndex: {
    appBar: 1200,
    drawer: 1100,
  },
  typography: {
    fontSize: 12,
  },
});
