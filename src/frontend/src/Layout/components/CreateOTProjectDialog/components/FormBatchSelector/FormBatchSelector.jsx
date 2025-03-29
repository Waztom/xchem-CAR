import { FormControl, FormHelperText, FormLabel, makeStyles } from '@material-ui/core';
import { ErrorMessage, useField } from 'formik';
import React from 'react';
import { BatchSelector } from '../../../../../common/components/BatchSelector';
import { SuspenseWithBoundary } from '../../../../../common/components/SuspenseWithBoundary';
import { OTWarningSection } from '../OTWarningSection';

const useStyles = makeStyles(theme => ({
  container: {
    display: 'grid',
    gap: theme.spacing()
  }
}));

export const FormBatchSelector = ({ name, label }) => {
  const classes = useStyles();

  const [field, meta, helpers] = useField(name);
  const selectedBatchesMap = field.value;

  return (
    <FormControl variant="filled" error={meta.touched && !!meta.error}>
      <FormLabel>{label}</FormLabel>

      <SuspenseWithBoundary>
        <div className={classes.container}>
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
        </div>
      </SuspenseWithBoundary>

      <ErrorMessage name={name}>{error => <FormHelperText error={true}>{error}</FormHelperText>}</ErrorMessage>
    </FormControl>
  );
};
