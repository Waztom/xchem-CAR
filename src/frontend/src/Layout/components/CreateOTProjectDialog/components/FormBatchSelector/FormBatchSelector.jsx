import React, { useEffect, useState } from 'react';
import { FormControl, FormLabel, FormHelperText, Box, Typography, Divider } from '@mui/material';
import { useField, useFormikContext } from 'formik';
import { BatchSelector } from '../../../../../common/components/BatchSelector';
import { StartingMaterialsUpload } from '../StartingMaterialsUpload/StartingMaterialsUpload';
import { DownloadTemplateButton } from '../DownloadTemplateButton/DownloadTemplateButton';

export const FormBatchSelector = ({ name, label }) => {
  const [field, meta, helpers] = useField(name);
  const [filesField] = useField('startingMaterialFiles'); // Get the files field
  const { setFieldValue } = useFormikContext();
  const [selectedBatches, setSelectedBatches] = useState([]);
  
  // Extract selected batch IDs whenever field.value changes
  useEffect(() => {
    const selected = Object.entries(field.value || {})
      .filter(([_, isSelected]) => isSelected)
      .map(([batchId]) => Number(batchId));
    
    setSelectedBatches(selected);
  }, [field.value]);

  return (
    <FormControl variant="filled" error={meta.touched && !!meta.error}>
      <FormLabel>{label}</FormLabel>

      <BatchSelector
        selectedBatchesMap={field.value}
        onBatchSelected={(batchId, selected) => {
          helpers.setValue({ 
            ...field.value, 
            [batchId]: selected 
          });
        }}
      />

      {selectedBatches.length > 0 && (
        <Box sx={{ mt: 3 }}>
          <Divider />
          <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mt: 2, mb: 1 }}>
            <Typography variant="subtitle2">
              Optional: Upload Starting Materials for Selected Batches
            </Typography>
            <DownloadTemplateButton />
          </Box>
          
          <Typography variant="caption" color="text.secondary">
            You can optionally provide custom starting material plate information for each batch
          </Typography>
          
          {selectedBatches.map(batchId => (
            <Box key={batchId} sx={{ mt: 2, mb: 2, pl: 2, borderLeft: '2px solid #e0e0e0' }}>
              <Typography variant="body2" fontWeight="medium">
                Batch ID: {batchId}
              </Typography>
              <StartingMaterialsUpload
                batchId={batchId}
                onFileChange={(batchId, file) => {
                  // Correctly update startingMaterialFiles with existing values
                  console.log(`Setting file for batch ${batchId}:`, file?.name || 'none');
                  setFieldValue('startingMaterialFiles', {
                    ...filesField.value, // Use the actual files field value
                    [batchId]: file
                  });
                }}
              />
            </Box>
          ))}
        </Box>
      )}

      {meta.touched && meta.error && (
        <FormHelperText error>{meta.error}</FormHelperText>
      )}
    </FormControl>
  );
};

FormBatchSelector.displayName = 'FormBatchSelector';
