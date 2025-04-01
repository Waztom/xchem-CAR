import React from 'react';
import { styled } from '@mui/material/styles';

const StyledIcon = styled('span')(({ theme }) => ({
  display: 'inline-flex',
  alignItems: 'center',
  justifyContent: 'center',
  width: 24,
  height: 24,
  '& > svg': {
    width: '100%',
    height: '100%',
    color: 'inherit'
  }
}));

export const IconComponent = ({ Component, ...rest }) => (
  <StyledIcon>
    <Component {...rest} />
  </StyledIcon>
);

IconComponent.displayName = 'IconComponent';
