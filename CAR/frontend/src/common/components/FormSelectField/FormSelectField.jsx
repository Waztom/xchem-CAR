import React from 'react';
import { ErrorMessage, Field, useField } from 'formik';
import { Select } from 'formik-material-ui';
import { FormControl, FormHelperText, InputLabel, MenuItem } from '@material-ui/core';

export const FormSelectField = ({ id, name, label, items, ...rest }) => {
  const meta = useField(name)[1];

  const labelId = `${id}-label`;

  return (
    <FormControl variant="outlined" error={meta.touched && !!meta.error}>
      <InputLabel id={labelId}>{label}</InputLabel>
      <Field
        component={Select}
        id={id}
        name={name}
        labelId={labelId}
        label={label}
        fullWidth
        {...rest}
        validate={() => false}
      >
        {items.map(({ value, text }) => (
          <MenuItem key={value} value={value}>
            {text}
          </MenuItem>
        ))}
      </Field>
      <ErrorMessage name={name}>{error => <FormHelperText error={true}>{error}</FormHelperText>}</ErrorMessage>
    </FormControl>
  );
};
