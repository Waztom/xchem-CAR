import { Button } from '@mui/material';
import { styled } from '@mui/material/styles';

const StyledButton = styled(Button)(({ theme }) => ({
  color: theme.palette.primary.contrastText
}));

export const SnackbarButton = props => (
  <StyledButton 
    variant="outlined" 
    color="inherit" 
    {...props} 
  />
);

SnackbarButton.displayName = 'SnackbarButton';
