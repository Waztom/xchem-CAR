import React from 'react';
import { styled } from '@mui/material/styles';
import { CircularProgress } from '@mui/material';
import { SuspenseWithBoundary } from '../SuspenseWithBoundary';

const StyledSection = styled('section')(({ theme }) => ({
  display: 'grid',
  gap: theme.spacing(2),
  padding: theme.spacing(2),
  position: 'relative',
  minHeight: 100
}));

const LoadingOverlay = styled('div')(({ theme }) => ({
  position: 'absolute',
  top: 0,
  left: 0,
  right: 0,
  bottom: 0,
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'center',
  backgroundColor: 'rgba(255, 255, 255, 0.7)',
  zIndex: theme.zIndex.modal - 1
}));

const LoadingFallback = () => (
  <LoadingOverlay>
    <CircularProgress />
  </LoadingOverlay>
);

const ErrorFallback = ({ error }) => (
  <div role="alert">
    <p>Something went wrong:</p>
    <pre style={{ color: 'red' }}>{error.message}</pre>
  </div>
);

export const DialogSection = ({ children, SuspenseProps, ErrorBoundaryProps }) => {
  return (
    <StyledSection>
      <SuspenseWithBoundary
        fallback={<LoadingFallback />}
        ErrorFallbackComponent={ErrorFallback}
        SuspenseProps={SuspenseProps}
        ErrorBoundaryProps={ErrorBoundaryProps}
      >
        {children}
      </SuspenseWithBoundary>
    </StyledSection>
  );
};

DialogSection.displayName = 'DialogSection';
