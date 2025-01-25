import { Switch, FormControlLabel } from '@material-ui/core';
import React from 'react';

export const LayoutSwitch = ({ checked, onChange, label }) => {
  return (
    <FormControlLabel
      control={<Switch checked={checked} onChange={(_, checked) => onChange(checked)} />}
      label={label}
    />
  );
};
