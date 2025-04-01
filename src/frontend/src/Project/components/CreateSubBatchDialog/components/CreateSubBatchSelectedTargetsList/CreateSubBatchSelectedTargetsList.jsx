import React, { memo } from 'react';
import { List, ListItem, Typography } from '@mui/material';
import { styled } from '@mui/material/styles';
import { useSelectedTargets } from './hooks/useSelectedTargets';
import { SuspenseWithBoundary } from '../../../../../common/components/SuspenseWithBoundary';

const StyledListItem = styled(ListItem)(({ theme }) => ({
  paddingTop: 0,
  paddingBottom: 0,
  gap: theme.spacing()
}));

const StyledImage = styled('img')({
  width: 120,
  height: 60,
  objectFit: 'contain'
});

/**
 * Renders a list of selected targets. Wrapped in a memo since it doesn't accept props and it's pointless to rerender
 * it each time a field in the Create subbatch dialog changes. Saves around 20ms on change in dev environment.
 */
const TargetsList = memo(() => {
  const selectedTargets = useSelectedTargets();

  return (
    <List disablePadding>
      {selectedTargets.map(({ target, methodsCount }) => (
        <StyledListItem key={target.id}>
          <StyledImage 
            src={target.image} 
            alt={target.name} 
          />
          <Typography variant="caption">
            <strong>{target.name}</strong>
            &nbsp;({methodsCount})
          </Typography>
        </StyledListItem>
      ))}
    </List>
  );
});

export const CreateSubBatchSelectedTargetsList = () => (
  <SuspenseWithBoundary>
    <TargetsList />
  </SuspenseWithBoundary>
);

CreateSubBatchSelectedTargetsList.displayName = 'CreateSubBatchSelectedTargetsList';
TargetsList.displayName = 'TargetsList';
