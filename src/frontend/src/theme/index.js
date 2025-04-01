import { createTheme, adaptV4Theme } from '@mui/material/styles';
import palette from './palette';

export const getTheme = () => {
  return createTheme(adaptV4Theme({
    palette,
    zIndex: {
      appBar: 1200,
      drawer: 1100
    },
    typography: {
      fontSize: 12
    }
  }));
};
