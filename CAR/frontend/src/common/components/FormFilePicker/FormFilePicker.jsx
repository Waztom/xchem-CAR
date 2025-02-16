import React from 'react';
import { ErrorMessage, useField } from 'formik';
import { Button, FormControl, FormHelperText, FormLabel, makeStyles } from '@material-ui/core';

const useStyles = makeStyles(theme => ({
  input: {
    display: 'none'
  }
}));

export const FormFilePicker = ({ name, label, description, id, buttonText, accept }) => {
  const classes = useStyles();

  const [field, meta, helpers] = useField(name);

  return (
    <FormControl variant="filled" error={meta.touched && !!meta.error}>
      <FormLabel>{label}</FormLabel>
      {description}
      <input
        accept={accept}
        className={classes.input}
        id={id}
        type="file"
        onChange={event => {
          const file = event.target.files[0];
          helpers.setValue(file || null);
          // Without timeout the setTouched would execute sooner than setValue
          setTimeout(() => {
            helpers.setTouched();
          });
        }}
        onBlur={() => helpers.setTouched()}
      />
      <label htmlFor={id}>
        <Button variant="contained" color="primary" component="span" fullWidth>
          {buttonText ?? 'Select file'}
        </Button>
      </label>
      <ErrorMessage name={name}>{error => <FormHelperText error={true}>{error}</FormHelperText>}</ErrorMessage>
      {!!field.value && <FormHelperText>{field.value?.name}</FormHelperText>}
    </FormControl>
  );
};
