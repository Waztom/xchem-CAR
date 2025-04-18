import { colors } from '@mui/material';

const white = '#FFFFFF';
const black = '#000000';
const darkBlack = '#000000de';

export default {
  darkBlack,
  black,
  white,
  primary: {
    contrastText: white,
    dark: colors.indigo[900],
    main: colors.indigo[500],
    light: colors.indigo[100],
    semidark: colors.indigo[300]
  },
  secondary: {
    contrastText: white,
    dark: colors.blue[900],
    main: colors.blue['A400'],
    light: colors.blue['A400']
  },
  success: {
    contrastText: white,
    dark: colors.green[900],
    main: colors.green[600],
    light: colors.green[400],
    lighter: colors.green[100]
  },
  info: {
    contrastText: white,
    dark: colors.blue[900],
    main: colors.blue[600],
    light: colors.blue[400]
  },
  warning: {
    contrastText: white,
    dark: colors.orange[900],
    darkLight: colors.orange[800],
    main: colors.orange[600],
    light: colors.orange[400]
  },
  error: {
    contrastText: white,
    dark: colors.red[900],
    main: colors.red[600],
    light: colors.red[400],
    lighter: colors.red[100]
  },
  text: {
    primary: colors.blueGrey[900],
    secondary: colors.blueGrey[600],
    link: colors.blue[600],
    disabled: colors.grey[500]
  },
  background: {
    default: '#F4F6F8',
    paper: white,
    divider: '#CECECE'
  },
  tag: {
    default: colors.grey[300]
  },
  icon: colors.blueGrey[600],
  divider: colors.grey[200],
  dividerDark: colors.grey[500]
};
