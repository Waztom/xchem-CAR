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

  // Required columns for the CSV file
  const requiredColumns = ["plate-ID", "labware-type", "well-index", "SMILES", "amount-uL"];

  // Validate CSV headers
  const validateCSVHeaders = (csvContent) => {
    try {
      // Get the first line and parse headers
      const firstLine = csvContent.split('\n')[0];
      const headers = firstLine.split(',').map(h => h.trim().toLowerCase());
      
      // Check if all required columns are present
      const missingColumns = requiredColumns.filter(col => 
        !headers.includes(col.toLowerCase())
      );
      
      if (missingColumns.length > 0) {
        return {
          valid: false,
          missingColumns
        };
      }
      
      return { valid: true };
    } catch (err) {
      return {
        valid: false,
        error: "Could not parse CSV headers"
      };
    }
  };

  const handleFileChange = (event) => {
    const file = event.target.files[0];
    if (!file) return;

    if (!file.name.endsWith('.csv')) {
      setError('File must be a CSV');
      setFileName(null);
      onFileChange(batchId, null);
      return;
    }

    // Read the file to validate headers
    const reader = new FileReader();
    reader.onload = (e) => {
      const content = e.target.result;
      const validation = validateCSVHeaders(content);
      
      if (!validation.valid) {
        if (validation.missingColumns) {
          setError(`Required columns missing: ${validation.missingColumns.join(', ')}`);
        } else {
          setError(validation.error || 'Invalid CSV format');
        }
        setFileName(null);
        onFileChange(batchId, null);
        return;
      }

      // CSV is valid
      setFileName(file.name);
      setError(null);
      onFileChange(batchId, file);
    };
    
    reader.onerror = () => {
      setError('Failed to read the file');
      setFileName(null);
      onFileChange(batchId, null);
    };
    
    reader.readAsText(file);
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
        <Box sx={{ mt: 1, mb: 0.5 }}>
          <Typography variant="body2" color="textSecondary" sx={{ mb: 1 }}>
            If provided, CSV must include: plate-ID, labware-type, well-index, SMILES, amount-uL
          </Typography>
        </Box>
      )}

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

      {fileName && (
        <Box sx={{ mt: 1 }}>
          <Typography variant="caption" color="textSecondary">
            <strong>Note:</strong> Each plate-ID/labware_type combination creates a separate physical plate
          </Typography>
        </Box>
      )}
    </Box>
  );
};