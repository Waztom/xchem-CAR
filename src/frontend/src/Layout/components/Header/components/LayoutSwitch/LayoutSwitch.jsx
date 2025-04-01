import { Switch, FormControlLabel } from '@mui/material';
import { styled } from '@mui/material/styles';
import React, { useState, useEffect } from 'react';

const ColorSwitch = styled(Switch)(({ theme, isChecked }) => ({
  '& .MuiSwitch-switchBase': {
    color: '#e91e63', // Default red-pink
    padding: 8, // Increased padding for larger hover area
    '&.Mui-checked': {
      color: '#4caf50', // Green when checked
      '&:hover': {
        backgroundColor: 'rgba(76, 175, 80, 0.15)' // Increased green opacity
      }
    },
    '&:hover': {
      backgroundColor: 'rgba(233, 30, 99, 0.15)' // Increased red-pink opacity
    }
  },
  '& .MuiSwitch-switchBase.Mui-checked + .MuiSwitch-track': {
    backgroundColor: '#4caf50', // Green track when checked
    opacity: 0.9 // Increased track opacity
  },
  '& .MuiSwitch-track': {
    backgroundColor: '#e91e63', // Red-pink track by default
    opacity: 0.9 // Increased track opacity
  }
}));

export const LayoutSwitch = ({ checked, onChange, label }) => {
  const [isChecked, setIsChecked] = useState(checked);

  useEffect(() => {
    setIsChecked(checked);
  }, [checked]);

  const handleChange = (event, checked) => {
    setIsChecked(checked);
    onChange(checked);
  };

  return (
    <FormControlLabel
      control={
        <ColorSwitch 
          checked={isChecked}
          onChange={handleChange}
        />
      }
      label={label}
    />
  );
};
