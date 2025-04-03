import React from 'react';
import { ErrorMessage, useField } from 'formik';
import { Button, FormControl, FormHelperText, FormLabel } from '@mui/material';
import { styled } from '@mui/material/styles';

const HiddenInput = styled('input')({
  display: 'none'
});

const FilePickerControl = styled(FormControl)(({ theme }) => ({
  display: 'flex',
  gap: theme.spacing(1)
}));

export const FormFilePicker = ({ 
  name, 
  label, 
  description, 
  id, 
  buttonText, 
  accept 
}) => {
  const [field, meta, helpers] = useField(name);

  return (
    <FilePickerControl variant="filled" error={meta.touched && !!meta.error}>
      <FormLabel>{label}</FormLabel>
      {description}
      <HiddenInput
        accept={accept}
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
        <Button 
          variant="contained" 
          color="primary" 
          component="span" 
          fullWidth
        >
          {buttonText ?? 'Select file'}
        </Button>
      </label>
      <ErrorMessage name={name}>
        {error => (
          <FormHelperText error>{error}</FormHelperText>
        )}
      </ErrorMessage>
      {!!field.value && (
        <FormHelperText>
          {field.value?.name}
        </FormHelperText>
      )}
    </FilePickerControl>
  );
};

FormFilePicker.displayName = 'FormFilePicker';
