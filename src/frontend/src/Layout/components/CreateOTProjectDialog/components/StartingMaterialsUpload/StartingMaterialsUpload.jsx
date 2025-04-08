import React, { useState } from 'react';
import { Button, Typography, Alert, Box, Tooltip } from '@mui/material';
import { CloudUpload, Delete } from '@mui/icons-material';
import { styled } from '@mui/material/styles';

const VisuallyHiddenInput = styled('input')({
  clip: 'rect(0 0 0 0)',
  clipPath: 'inset(50%)',
  height: 1,
  overflow: 'hidden',
  position: 'absolute',
  bottom: 0,
  left: 0,
  whiteSpace: 'nowrap',
  width: 1,
});

export const StartingMaterialsUpload = ({ batchId, onFileChange }) => {
  const [error, setError] = useState(null);
  const [fileName, setFileName] = useState(null);

  const handleFileChange = (event) => {
    const file = event.target.files[0];
    if (!file) return;

    if (!file.name.endsWith('.csv')) {
      setError('File must be a CSV');
      setFileName(null);
      onFileChange(batchId, null);
      return;
    }

    setFileName(file.name);
    setError(null);
    onFileChange(batchId, file);
  };

  const handleRemoveFile = () => {
    setFileName(null);
    setError(null);
    onFileChange(batchId, null);
  };

  return (
    <Box sx={{ mt: 1, mb: 1 }}>
      <Tooltip title="Upload custom starting material plate information (optional)">
        <Button
          component="label"
          variant="outlined"
          startIcon={<CloudUpload />}
          size="small"
        >
          Upload CSV (Optional)
          <VisuallyHiddenInput 
            type="file" 
            accept=".csv"
            onChange={handleFileChange}
          />
        </Button>
      </Tooltip>

      {fileName && (
        <Box sx={{ display: 'flex', alignItems: 'center', mt: 0.5 }}>
          <Typography variant="caption" sx={{ flexGrow: 1 }}>
            Selected: {fileName}
          </Typography>
          <Button 
            size="small" 
            color="error" 
            startIcon={<Delete />}
            onClick={handleRemoveFile}
          >
            Remove
          </Button>
        </Box>
      )}

      {error && (
        <Alert severity="error" sx={{ mt: 1 }}>
          {error}
        </Alert>
      )}
    </Box>
  );
};