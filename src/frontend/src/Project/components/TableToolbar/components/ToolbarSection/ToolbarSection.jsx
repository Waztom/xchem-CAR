import React from 'react';
import { Typography } from '@mui/material';
import { styled } from '@mui/material/styles';

const StyledSection = styled('section')(({ theme }) => ({
  marginBottom: theme.spacing()
}));

const ContentWrapper = styled('div')(({ theme }) => ({
  display: 'grid',
  gap: theme.spacing()
}));

export const ToolbarSection = ({ title, children }) => (
  <StyledSection>
    <Typography variant="h6" component="p">
      {title}
    </Typography>
    <ContentWrapper>{children}</ContentWrapper>
  </StyledSection>
);

ToolbarSection.displayName = 'ToolbarSection';
