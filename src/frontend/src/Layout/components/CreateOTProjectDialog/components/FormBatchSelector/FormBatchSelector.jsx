import { FormControl, FormHelperText, FormLabel } from '@mui/material';
import { styled } from '@mui/material/styles';
import { ErrorMessage, useField } from 'formik';
import React from 'react';
import { BatchSelector } from '../../../../../common/components/BatchSelector';
import { SuspenseWithBoundary } from '../../../../../common/components/SuspenseWithBoundary';
import { OTWarningSection } from '../OTWarningSection';

const StyledContainer = styled('div')(({ theme }) => ({
  display: 'grid',
  gap: theme.spacing()
}));

const StyledFormControl = styled(FormControl)(({ theme }) => ({
  '& .MuiFormLabel-root': {
    marginBottom: theme.spacing(1)
  }
}));

export const FormBatchSelector = ({ name, label }) => {
  const [field, meta, helpers] = useField(name);
  const selectedBatchesMap = field.value;

  return (
    <StyledFormControl variant="filled" error={meta.touched && !!meta.error}>
      <FormLabel>{label}</FormLabel>

      <SuspenseWithBoundary>
        <StyledContainer>
          <BatchSelector
            selectedBatchesMap={selectedBatchesMap}
            onBatchSelected={(batchId, selected) => {
              helpers.setValue({ ...selectedBatchesMap, [batchId]: selected });
              // Without timeout the setTouched would execute sooner than setValue
              setTimeout(() => {
                helpers.setTouched();
              });
            }}
          />
          <OTWarningSection selectedBatchesMap={selectedBatchesMap} />
        </StyledContainer>
      </SuspenseWithBoundary>

      <ErrorMessage name={name}>
        {error => <FormHelperText error>{error}</FormHelperText>}
      </ErrorMessage>
    </StyledFormControl>
  );
};

FormBatchSelector.displayName = 'FormBatchSelector';
