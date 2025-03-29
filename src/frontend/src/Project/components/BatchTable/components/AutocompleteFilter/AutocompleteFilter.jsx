import React from 'react';
import { Autocomplete } from '@material-ui/lab';
import { Chip, makeStyles, TextField } from '@material-ui/core';

const useStyles = makeStyles(theme => ({
  chip: {
    maxWidth: 240
  }
}));

export const AutocompleteFilter = ({ id, options, label, placeholder, filterValue, setFilter }) => {
  const classes = useStyles();

  return (
    <Autocomplete
      multiple
      fullWidth
      id={id}
      options={options}
      getOptionLabel={option => String(option)}
      renderOption={option => (!!option ? option : <i>{String(option)}</i>)}
      filterSelectedOptions
      renderInput={params => <TextField {...params} variant="outlined" label={label} placeholder={placeholder} />}
      renderTags={(values, getTagProps) =>
        values.map((value, index) => (
          <Chip
            {...getTagProps({ index })}
            classes={{ root: classes.chip }}
            label={!!value ? value : <i>{String(value)}</i>}
          />
        ))
      }
      value={filterValue || []}
      onChange={(_, value) => setFilter(value)}
    />
  );
};
