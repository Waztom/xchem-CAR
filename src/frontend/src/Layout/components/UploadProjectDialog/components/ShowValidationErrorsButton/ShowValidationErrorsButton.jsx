import { useSnackbar } from 'notistack';
import React from 'react';
import { SnackbarButton } from '../../../../../common/components/SnackbarButton';
import { requestShowSmilesValidationErrors } from '../../../../stores/smilesValidationErrorsDialogStore';

export const ShowValidationErrorsButton = ({ messageId, errors }) => {
  const { closeSnackbar } = useSnackbar();

  return (
    <SnackbarButton
      onClick={() => {
        closeSnackbar(messageId);

        requestShowSmilesValidationErrors(errors);
      }}
    >
      Show details
    </SnackbarButton>
  );
};
