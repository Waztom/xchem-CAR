import React from 'react';
import { saveAs } from 'file-saver';
import { SnackbarButton } from '../../../../../common/components/SnackbarButton';
import { useSnackbar } from 'notistack';
import { useCurrentProjectStore } from '../../../../../common/stores/currentProjectStore';

export const DownloadButton = ({ messageId, csvData }) => {
  const { closeSnackbar } = useSnackbar();

  const project = useCurrentProjectStore.useCurrentProject();

  return (
    <SnackbarButton
      onClick={() => {
        const blob = new Blob([csvData], { type: 'text/plain' });
        saveAs(blob, `${project.name}-${new Date().toISOString()}.csv`);
        closeSnackbar(messageId);
      }}
    >
      Download
    </SnackbarButton>
  );
};
