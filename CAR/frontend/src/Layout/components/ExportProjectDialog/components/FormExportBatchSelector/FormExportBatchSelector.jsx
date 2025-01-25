import { FormControl, FormHelperText, FormLabel } from '@material-ui/core';
import { ErrorMessage, useField } from 'formik';
import React from 'react';
import { BatchSelector } from '../../../../../common/components/BatchSelector';
import { SuspenseWithBoundary } from '../../../../../common/components/SuspenseWithBoundary';

export const FormExportBatchSelector = ({ name, label }) => {
  const [field, meta, helpers] = useField(name);
  const selectedBatchesMap = field.value;

  return (
    <FormControl variant="filled" error={meta.touched && !!meta.error}>
      <FormLabel>{label}</FormLabel>

      <SuspenseWithBoundary>
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
      </SuspenseWithBoundary>

      <ErrorMessage name={name}>{error => <FormHelperText error={true}>{error}</FormHelperText>}</ErrorMessage>
    </FormControl>
  );
};
