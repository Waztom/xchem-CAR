import React from 'react';
import { List, ListItem, Typography, Alert } from '@mui/material';
import { styled } from '@mui/material/styles';
import { useGetIncompatibleTargets } from './hooks/useGetIncompatibleTargets';
import { SuspenseWithBoundary } from '../../../../../common/components/SuspenseWithBoundary';

const StyledListItem = styled(ListItem)(({ theme }) => ({
  paddingTop: 0,
  paddingBottom: 0,
  gap: theme.spacing()
}));

const StyledImage = styled('img')({
  mixBlendMode: 'multiply'
});

const WarningContent = ({ selectedBatchesMap }) => {
  const incompatibleTargets = useGetIncompatibleTargets(selectedBatchesMap);

  return (
    <>
      {incompatibleTargets.map(({ batch, targets }) => (
        <Alert key={batch.id} severity="warning" variant="outlined">
          These targets from batch <strong>{batch.batchtag}</strong> contain only methods which can&apos;t be executed on
          OpenTrons:
          <List disablePadding>
            {targets.map(target => (
              <StyledListItem key={target.id}>
                <StyledImage 
                  src={target.image} 
                  width={120} 
                  height={60} 
                  alt={target.name} 
                />
                <Typography variant="caption">
                  <strong>{target.name}</strong>
                </Typography>
              </StyledListItem>
            ))}
          </List>
        </Alert>
      ))}
    </>
  );
};

export const OTWarningSection = (props) => (
  <SuspenseWithBoundary>
    <WarningContent {...props} />
  </SuspenseWithBoundary>
);

OTWarningSection.displayName = 'OTWarningSection';
