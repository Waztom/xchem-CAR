import React, { useState } from 'react';
import { Autocomplete, TextField, Chip } from '@mui/material';
import { styled } from '@mui/material/styles';

const StyledChip = styled(Chip)(({ theme }) => ({
  maxWidth: 240,
  '& .MuiChip-label': {
    whiteSpace: 'nowrap',
    overflow: 'hidden',
    textOverflow: 'ellipsis'
  }
}));

export const AutocompleteFilter = ({ 
  id, 
  options, 
  label, 
  placeholder, 
  filterValue, 
  setFilter
}) => {
  const [inputValue, setInputValue] = useState('');

  return (
    <Autocomplete
      multiple
      fullWidth
      freeSolo
      id={id}
      options={options}
      inputValue={inputValue}
      onInputChange={(event, newValue) => {
        setInputValue(newValue);
      }}
      getOptionLabel={option => String(option)}
      renderOption={(props, option) => (
        <li {...props}>
          {!!option ? option : <i>{String(option)}</i>}
        </li>
      )}
      filterSelectedOptions
      renderInput={params => (
        <TextField 
          {...params} 
          variant="outlined" 
          label={label} 
          placeholder={placeholder}
        />
      )}
      renderTags={(values, getTagProps) =>
        values.map((value, index) => (
          <StyledChip
            key={value}
            {...getTagProps({ index })}
            label={!!value ? value : <i>{String(value)}</i>}
          />
        ))
      }
      value={filterValue || []}
      onChange={(_, value) => {
        setFilter(value);
      }}
    />
  );
};

AutocompleteFilter.displayName = 'AutocompleteFilter';
