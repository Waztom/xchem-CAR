import React from 'react';
import { useSnackbar } from 'notistack';
import { requestOtProjectSummary } from '../../../../stores/otProjectSummaryDialogStore';
import { SnackbarButton } from '../../../../../common/components/SnackbarButton';

export const ShowOTProjectSummaryButton = ({ messageId, otProjectId }) => {
  const { closeSnackbar } = useSnackbar();

  return (
    <SnackbarButton
      onClick={() => {
        closeSnackbar(messageId);
        requestOtProjectSummary(otProjectId);
      }}
    >
      Show summary
    </SnackbarButton>
  );
};
