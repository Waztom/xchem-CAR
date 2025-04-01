import React, { useLayoutEffect, useState } from 'react';
import { FormControl, FormLabel, IconButton, Slider, Tooltip } from '@mui/material';
import { styled } from '@mui/material/styles';
import { Cancel } from '@mui/icons-material';

const FilterWrapper = styled('div')(({ theme }) => ({
  marginBottom: theme.spacing(2),
  display: 'flex',
  alignItems: 'flex-start'
}));

const StyledSlider = styled(Slider)(({ theme }) => ({
  margin: `0 ${theme.spacing()}px`,
  width: `calc(100% - ${theme.spacing(2)}px)`
}));

const ClearButtonWrapper = styled('div')(({ theme }) => ({
  paddingTop: theme.spacing(0.625)
}));

export const RangeFilter = ({ id, label, min, max, filterValue, setFilter }) => {
  const labelId = `${id}-label`;
  const [value, setValue] = useState(filterValue);

  useLayoutEffect(() => {
    setValue(filterValue);
  }, [filterValue]);

  // In case there are no numbers to take min and max from the STD returns +/-Infinity
  const valid = Number.isFinite(min) && Number.isFinite(max);
  const active = !!filterValue;
  const clearEnabled = valid && active;

  return (
    <FilterWrapper>
      <FormControl fullWidth>
        <FormLabel>{label}</FormLabel>
        <StyledSlider
          value={value || [min, max]}
          onChange={(_, value) => setValue(Array.isArray(value) ? value : [value, value])}
          aria-labelledby={labelId}
          id={id}
          getAriaValueText={value => `${label} - ${value}`}
          min={min}
          max={max}
          marks={
            valid
              ? [min, max, ...(value || [])]
                  .reduce((marks, nextValue) => {
                    if (!marks.includes(nextValue)) {
                      marks.push(nextValue);
                    }
                    return marks;
                  }, [])
                  .map(val => ({ value: val, label: val }))
              : []
          }
          onChangeCommitted={(_, value) => setFilter(value)}
          disabled={!valid}
          color={active ? 'primary' : 'secondary'}
        />
      </FormControl>
      <ClearButtonWrapper>
        <Tooltip title={clearEnabled ? 'Clear filter' : ''}>
          <span>
            <IconButton 
              onClick={() => setFilter()} 
              disabled={!clearEnabled} 
              size="large"
            >
              <Cancel />
            </IconButton>
          </span>
        </Tooltip>
      </ClearButtonWrapper>
    </FilterWrapper>
  );
};

RangeFilter.displayName = 'RangeFilter';
