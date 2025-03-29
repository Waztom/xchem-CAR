import React from 'react';
import { Chip, FormControl, FormHelperText, IconButton, makeStyles, TextField, Tooltip } from '@material-ui/core';
import { ErrorMessage, useField } from 'formik';
import { Autocomplete } from '@material-ui/lab';
import { createFilterOptions } from '@material-ui/lab/Autocomplete';
import classNames from 'classnames';
import { CloudUpload } from '@material-ui/icons';

const useStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
    alignItems: 'center',
    gap: theme.spacing()
  },
  chip: {
    maxWidth: 240
  },
  error: {
    backgroundColor: theme.palette.error.light,
    color: theme.palette.error.contrastText
  },
  errorIcon: {
    color: theme.palette.error.dark
  },
  input: {
    display: 'none'
  }
}));

export const FormReactionIdsSelector = ({ name, label, autocompleteId, filePickerId, reactionsIds, reactionsMap }) => {
  const classes = useStyles();

  const [field, meta, helpers] = useField(name);

  return (
    <FormControl variant="filled" error={meta.touched && !!meta.error} fullWidth>
      <div className={classes.root}>
        <Autocomplete
          multiple
          fullWidth
          id={autocompleteId}
          options={reactionsIds}
          getOptionLabel={option => String(option.id)}
          renderOption={option => (
            <>
              {option}&nbsp;<i>{reactionsMap[option].reactionclass}</i>
            </>
          )}
          filterOptions={createFilterOptions({
            stringify: option => `${option} ${reactionsMap[option].reactionclass}`
          })}
          filterSelectedOptions
          renderInput={params => <TextField {...params} variant="outlined" label={label} placeholder={label} />}
          renderTags={(values, getTagProps) =>
            values.map((value, index) => {
              return (
                <Chip
                  {...getTagProps({ index })}
                  classes={{
                    root: classNames(classes.chip, !reactionsMap[value] && classes.error),
                    deleteIcon: classNames(!reactionsMap[value] && classes.errorIcon)
                  }}
                  label={value}
                />
              );
            })
          }
          value={field.value}
          onChange={(_, value) => helpers.setValue(value)}
        />
        <input
          accept="text/csv"
          className={classes.input}
          id={filePickerId}
          type="file"
          onChange={event => {
            const file = event.target.files[0];

            if (file) {
              const reader = new FileReader();
              reader.addEventListener('load', () => {
                const ids = reader.result.split(';');
                helpers.setValue([...new Set([...field.value, ...ids])]);
              });

              reader.readAsText(file);
            }
          }}
        />
        <label htmlFor={filePickerId}>
          <Tooltip title="Upload reaction IDs from file">
            <IconButton component="span">
              <CloudUpload />
            </IconButton>
          </Tooltip>
        </label>
      </div>
      <ErrorMessage name={name}>{error => <FormHelperText error={true}>{error}</FormHelperText>}</ErrorMessage>
    </FormControl>
  );
};
