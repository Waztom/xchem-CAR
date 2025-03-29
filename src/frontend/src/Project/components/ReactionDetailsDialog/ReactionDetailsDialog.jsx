import { makeStyles, Typography } from '@material-ui/core';
import React from 'react';
import { CloseDialog } from '../../../common/components/CloseDialog';
import { DialogSection } from '../../../common/components/DialogSection';
import { DialogSectionHeading } from '../../../common/components/DialogSectionHeading';
import {
  setReactionDetailsDialogOpen,
  setReactionForReactionDetailsDialog,
  useReactionDetailsDialogStore
} from '../../stores/reactionDetailsDialogStore';
import { ProductSection } from './components/ProductSection/ProductSection';
import { ReactantSection } from './components/ReactantSection/ReactantSection';

const useStyles = makeStyles(theme => ({
  root: {
    display: 'grid',
    gap: theme.spacing(2),
    '& > :nth-child(even)': {
      backgroundColor: theme.palette.action.hover
    }
  },
  reactionWrapper: {
    display: 'grid',
    justifyContent: 'center'
  },
  reactionNameWrapper: {
    textAlign: 'center',
    display: 'grid'
  },
  image: {
    mixBlendMode: 'multiply'
  }
}));

export const ReactionDetailsDialog = () => {
  const classes = useStyles();

  const { dialogOpen, reaction } = useReactionDetailsDialogStore();

  return (
    <CloseDialog
      id="reaction-details-dialog"
      open={dialogOpen}
      title="Reaction details"
      maxWidth="md"
      fullWidth
      content={
        <div className={classes.root}>
          <DialogSection>
            <DialogSectionHeading>Reaction</DialogSectionHeading>
            <div className={classes.reactionWrapper}>
              <img className={classes.image} src={reaction?.image} width={270} height={60} />
              <div className={classes.reactionNameWrapper}>
                <Typography noWrap>{reaction?.reactionclass}</Typography>
              </div>
            </div>
          </DialogSection>

          {reaction?.reactants.map((reactant, index) => (
            <ReactantSection key={reactant.id} reactant={reactant} index={index} />
          ))}

          {!!reaction && <ProductSection product={reaction?.products[0]} />}
        </div>
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
