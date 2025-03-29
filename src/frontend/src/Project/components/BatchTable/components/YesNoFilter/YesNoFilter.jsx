import { FormControl, InputLabel, MenuItem, Select } from '@material-ui/core';
import React from 'react';

export const YesNoFilter = ({ id, label, filterValue, setFilter }) => {
  const labelId = `${id}-label`;

  return (
    <FormControl variant="outlined" fullWidth>
      <InputLabel id={labelId}>{label}</InputLabel>
      <Select
        labelId={labelId}
        id={id}
        value={filterValue ?? ''}
        onChange={event => setFilter(event.target.value)}
        label={label}
      >
        <MenuItem value={''}>
          <i>None</i>
        </MenuItem>
        <MenuItem value={true}>Yes</MenuItem>
        <MenuItem value={false}>No</MenuItem>
      </Select>
    </FormControl>
  );
};
