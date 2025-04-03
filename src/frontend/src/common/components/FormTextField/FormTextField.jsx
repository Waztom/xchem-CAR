import { Field } from 'formik';
import { TextField } from 'formik-mui'; 
import React from 'react';

export const FormTextField = ({ name, label, placeholder, ...rest }) => {
  return (
    <Field
      component={TextField}
      label={label}
      name={name}
      placeholder={placeholder}
      variant="outlined"
      fullWidth
      {...rest}
    />
  );
};
