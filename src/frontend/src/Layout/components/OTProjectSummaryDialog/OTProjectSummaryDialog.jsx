import { Typography } from '@material-ui/core';
import React from 'react';
import { DialogSection } from '../../../common/components/DialogSection';
import { DialogSectionHeading } from '../../../common/components/DialogSectionHeading';
import { SuspenseWithBoundary } from '../../../common/components/SuspenseWithBoundary';
import { CloseDialog } from '../../../common/components/CloseDialog';
import {
  setOtProjectForSummaryDialog,
  setOtProjectSummaryDialogOpen,
  useOtProjectSummaryDialogStore
} from '../../stores/otProjectSummaryDialogStore';
import { BatchProtocolList } from './components/BatchProtocolList';

export const OTProjectSummaryDialog = () => {
  const { dialogOpen, otProjectId } = useOtProjectSummaryDialogStore();

  return (
    <CloseDialog
      id="ot-protocol-summary-dialog"
      open={dialogOpen}
      title="OT protocol summary"
      content={
        <DialogSection>
          <DialogSectionHeading>OT protocols</DialogSectionHeading>
          <Typography>
            This is a list of batches for which OT protocols have been generated with download links:
          </Typography>
          <SuspenseWithBoundary>
            <BatchProtocolList otProjectId={otProjectId} />
          </SuspenseWithBoundary>
        </DialogSection>
      }
      onClose={() => {
        setOtProjectSummaryDialogOpen(false);
      }}
      TransitionProps={{
        onExited: () => setOtProjectForSummaryDialog(null)
      }}
    />
  );
};
