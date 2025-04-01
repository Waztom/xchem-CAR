import React, { forwardRef } from 'react';
import { Paper, Typography } from '@mui/material';
import { styled } from '@mui/material/styles';
import { SuspenseWithBoundary } from '../SuspenseWithBoundary';

const TitleWrapper = styled('div')(({ theme }) => ({
  height: theme.spacing(6),
  padding: `0 ${theme.spacing(2)}`,
  display: 'flex',
  alignItems: 'center',
  color: theme.palette.common.white,
  backgroundColor: theme.palette.primary.main,
  gap: theme.spacing()
}));

const Title = styled(Typography)(() => ({
  flexGrow: 1
}));

export const ContentBox = forwardRef(
  ({ title, children, endAdornment, PaperProps = {}, SuspenseProps, ErrorBoundaryProps }, ref) => {
    return (
      <Paper ref={ref} square {...PaperProps}>
        <TitleWrapper>
          <Title variant="h6" component="h2" noWrap>
            {title}
          </Title>
          {endAdornment}
        </TitleWrapper>
        <SuspenseWithBoundary SuspenseProps={SuspenseProps} ErrorBoundaryProps={ErrorBoundaryProps}>
          {children}
        </SuspenseWithBoundary>
      </Paper>
    );
  }
);

ContentBox.displayName = 'ContentBox';
