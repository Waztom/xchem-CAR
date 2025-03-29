import { FormControl, FormLabel, IconButton, makeStyles, Slider, Tooltip } from '@material-ui/core';
import { Cancel } from '@material-ui/icons';
import React, { useLayoutEffect, useState } from 'react';

const useStyles = makeStyles(theme => ({
  root: {
    marginBottom: theme.spacing(2),
    display: 'flex',
    alignItems: 'flex-start'
  },
  slider: {
    margin: `0 ${theme.spacing()}px`,
    width: `calc(100% - ${theme.spacing(2)}px)`
  },
  clearWrapper: {
    paddingTop: 5
  }
}));

export const RangeFilter = ({ id, label, min, max, filterValue, setFilter }) => {
  const classes = useStyles();

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
    <div className={classes.root}>
      <FormControl fullWidth>
        <FormLabel>{label}</FormLabel>
        <Slider
          className={classes.slider}
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
      <div className={classes.clearWrapper}>
        <Tooltip title={clearEnabled ? 'Clear filter' : ''}>
          <IconButton onClick={() => setFilter()} disabled={!clearEnabled}>
            <Cancel />
          </IconButton>
        </Tooltip>
      </div>
    </div>
  );
};
