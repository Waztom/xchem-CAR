import React from 'react';
import { CircularProgress } from '@mui/material';
import { styled } from '@mui/material/styles';
import { SuspenseWithBoundary } from '../SuspenseWithBoundary';

const SpinnerWrapper = styled('div')(({ theme }) => ({
  width: '100%',
  display: 'grid',
  placeContent: 'center',
  padding: theme.spacing(2)
}));

const LoadingSpinnerContent = () => (
  <SpinnerWrapper>
    <CircularProgress />
  </SpinnerWrapper>
);

const ErrorFallback = ({ error }) => (
  <SpinnerWrapper>
    <div role="alert">
      <p>Error loading content:</p>
      <pre style={{ color: theme.palette.error.main }}>{error.message}</pre>
    </div>
  </SpinnerWrapper>
);

export const LoadingSpinner = () => (
  <SuspenseWithBoundary
    fallback={<CircularProgress />}
    ErrorFallbackComponent={ErrorFallback}
  >
    <LoadingSpinnerContent />
  </SuspenseWithBoundary>
);

LoadingSpinner.displayName = 'LoadingSpinner';
