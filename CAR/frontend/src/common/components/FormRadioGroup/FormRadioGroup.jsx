import React from 'react';
import { FormControl, FormControlLabel, FormLabel, FormHelperText, Radio } from '@material-ui/core';
import { ErrorMessage, Field, useField } from 'formik';
import { RadioGroup } from 'formik-material-ui';

export const FormRadioGroup = ({ name, label, options, ...rest }) => {
  const meta = useField(name)[1];

  return (
    <FormControl component="fieldset" variant="filled" error={meta.touched && !!meta.error}>
      <FormLabel component="legend">{label}</FormLabel>
      <Field component={RadioGroup} name={name} {...rest}>
        {options.map(({ value, label }) => (
          <FormControlLabel key={value} value={value} control={<Radio />} label={label} />
        ))}
      </Field>
      <ErrorMessage name={name}>{msg => <FormHelperText error={true}>{msg}</FormHelperText>}</ErrorMessage>
    </FormControl>
  );
};
