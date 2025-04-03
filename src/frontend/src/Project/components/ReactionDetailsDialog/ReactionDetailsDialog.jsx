import React from 'react';
import { Typography } from '@mui/material';
import { styled } from '@mui/material/styles';
import { CloseDialog } from '../../../common/components/CloseDialog';
import { DialogSection } from '../../../common/components/DialogSection';
import { DialogSectionHeading } from '../../../common/components/DialogSectionHeading';
import { SuspenseWithBoundary } from '../../../common/components/SuspenseWithBoundary';
import {
  setReactionDetailsDialogOpen,
  setReactionForReactionDetailsDialog,
  useReactionDetailsDialogStore
} from '../../stores/reactionDetailsDialogStore';
import { ProductSection } from './components/ProductSection/ProductSection';
import { ReactantSection } from './components/ReactantSection/ReactantSection';

const StyledRoot = styled('div')(({ theme }) => ({
  display: 'grid',
  gap: theme.spacing(2),
  '& > :nth-child(even)': {
    backgroundColor: theme.palette.action.hover
  }
}));

const ReactionWrapper = styled('div')(() => ({
  display: 'grid',
  justifyContent: 'center'
}));

const ReactionNameWrapper = styled('div')(() => ({
  textAlign: 'center',
  display: 'grid'
}));

const StyledImage = styled('img')(() => ({
  mixBlendMode: 'multiply'
}));

const ReactionDetailsContent = () => {
  const { reaction } = useReactionDetailsDialogStore();

  return (
    <StyledRoot>
      <DialogSection>
        <DialogSectionHeading>Reaction</DialogSectionHeading>
        <ReactionWrapper>
          <StyledImage src={reaction?.image} width={270} height={60} />
          <ReactionNameWrapper>
            <Typography noWrap>{reaction?.reactionclass}</Typography>
          </ReactionNameWrapper>
        </ReactionWrapper>
      </DialogSection>

      {reaction?.reactants.map((reactant, index) => (
        <ReactantSection key={reactant.id} reactant={reactant} index={index} />
      ))}

      {!!reaction && <ProductSection product={reaction?.products[0]} />}
    </StyledRoot>
  );
};

const LoadingFallback = () => (
  <StyledRoot>
    <DialogSection>
      <Typography>Loading reaction details...</Typography>
    </DialogSection>
  </StyledRoot>
);

export const ReactionDetailsDialog = () => {
  const { dialogOpen } = useReactionDetailsDialogStore();

  return (
    <CloseDialog
      id="reaction-details-dialog"
      open={dialogOpen}
      title="Reaction details"
      maxWidth="md"
      fullWidth
      content={
        <SuspenseWithBoundary fallback={<LoadingFallback />}>
          <ReactionDetailsContent />
        </SuspenseWithBoundary>
      }
      onClose={() => {
        setReactionDetailsDialogOpen(false);
      }}
      TransitionProps={{
        onExited: () => setReactionForReactionDetailsDialog(null)
      }}
    />
  );
};

ReactionDetailsDialog.displayName = 'ReactionDetailsDialog';
