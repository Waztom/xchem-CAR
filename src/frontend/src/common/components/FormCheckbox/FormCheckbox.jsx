import { Field } from 'formik';
import { CheckboxWithLabel } from 'formik-mui';
import React from 'react';

export const FormCheckbox = ({ name, label, ...rest }) => {
  return (
    <Field
      component={CheckboxWithLabel}
      type="checkbox"
      name={name}
      Label={{ label: label }}
      {...rest}
    />
  );
};
