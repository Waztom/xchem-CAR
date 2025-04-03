import React from 'react';
import { 
  Chip, 
  FormControl, 
  FormHelperText, 
  IconButton, 
  TextField, 
  Tooltip,
  Autocomplete,
  createFilterOptions 
} from '@mui/material';
import { styled } from '@mui/material/styles';
import { ErrorMessage, useField } from 'formik';
import { CloudUpload } from '@mui/icons-material';

const Container = styled('div')(({ theme }) => ({
  display: 'flex',
  alignItems: 'center',
  gap: theme.spacing()
}));

const StyledChip = styled(Chip, {
  shouldForwardProp: prop => prop !== 'isError'
})(({ theme, isError }) => ({
  maxWidth: 240,
  ...(isError && {
    backgroundColor: theme.palette.error.light,
    color: theme.palette.error.contrastText,
    '& .MuiChip-deleteIcon': {
      color: theme.palette.error.dark
    }
  })
}));

const HiddenInput = styled('input')({
  display: 'none'
});

export const FormReactionIdsSelector = ({ 
  name, 
  label, 
  autocompleteId, 
  filePickerId, 
  reactionsIds, 
  reactionsMap 
}) => {
  const [field, meta, helpers] = useField(name);

  const handleFileChange = (event) => {
    const file = event.target.files[0];
    if (file) {
      const reader = new FileReader();
      reader.addEventListener('load', () => {
        const ids = reader.result.split(';');
        helpers.setValue([...new Set([...field.value, ...ids])]);
      });
      reader.readAsText(file);
    }
  };

  return (
    <FormControl variant="filled" error={meta.touched && !!meta.error} fullWidth>
      <Container>
        <Autocomplete
          multiple
          fullWidth
          id={autocompleteId}
          options={reactionsIds}
          getOptionLabel={option => String(option.id)}
          renderOption={(props, option) => (
            <li {...props}>
              {option}&nbsp;<i>{reactionsMap[option].reactionclass}</i>
            </li>
          )}
          filterOptions={createFilterOptions({
            stringify: option => `${option} ${reactionsMap[option].reactionclass}`
          })}
          filterSelectedOptions
          renderInput={params => (
            <TextField 
              {...params} 
              variant="outlined" 
              label={label} 
              placeholder={label} 
            />
          )}
          renderTags={(values, getTagProps) =>
            values.map((value, index) => (
              <StyledChip
                {...getTagProps({ index })}
                key={value}
                label={value}
                isError={!reactionsMap[value]}
              />
            ))
          }
          value={field.value}
          onChange={(_, value) => helpers.setValue(value)}
        />
        <HiddenInput
          accept="text/csv"
          id={filePickerId}
          type="file"
          onChange={handleFileChange}
        />
        <label htmlFor={filePickerId}>
          <Tooltip title="Upload reaction IDs from file">
            <IconButton component="span" size="large">
              <CloudUpload />
            </IconButton>
          </Tooltip>
        </label>
      </Container>
      <ErrorMessage name={name}>
        {error => <FormHelperText error>{error}</FormHelperText>}
      </ErrorMessage>
    </FormControl>
  );
};

FormReactionIdsSelector.displayName = 'FormReactionIdsSelector';
