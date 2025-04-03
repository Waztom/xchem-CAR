import React from 'react';
import { Typography } from '@mui/material';
import { styled } from '@mui/material/styles';

const StyledTypography = styled(Typography)(({ theme }) => ({
  fontSize: '0.9rem',
  fontWeight: 500,
  marginBottom: theme.spacing(1)
}));

export const DialogSectionHeading = ({ children }) => (
  <StyledTypography variant="h6" component="h3">
    {children}
  </StyledTypography>
);

DialogSectionHeading.displayName = 'DialogSectionHeading';
