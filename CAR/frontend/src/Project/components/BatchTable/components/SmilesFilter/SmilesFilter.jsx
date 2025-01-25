import { Chip, IconButton, makeStyles, TextField, Tooltip } from '@material-ui/core';
import { AddCircle, Clear, CloudUpload } from '@material-ui/icons';
import classNames from 'classnames';
import React, { useState } from 'react';
import { CanonicalizeSmilesDialog } from '../../../CanonicalizeSmilesDialog';
import { useCanonicalizeSmiles } from '../../../../hooks/useCanonicalizeSmiles';

const useStyles = makeStyles(theme => ({
  root: {
    display: 'flex',
    alignItems: 'center',
    gap: theme.spacing()
  },
  input: {
    padding: 9,
    paddingRight: 65,
    flexWrap: 'wrap',
    '&:hover .smiles-filter-dirty': {
      visibility: 'visible'
    }
  },
  inputField: {
    flexGrow: 1,
    minWidth: 30,
    width: 0,
    padding: '9.5px 4px'
  },
  chip: {
    margin: 3,
    maxWidth: 240
  },
  button: {
    padding: 2,
    marginRight: -2
  },
  clear: {
    visibility: 'hidden'
  },
  clearActive: {
    visibility: 'visible'
  },
  endAdornment: {
    position: 'absolute',
    right: 9,
    top: 'calc(50% - 14px)'
  }
}));

export const SmilesFilter = ({ id, label, filterValue = [], setFilter }) => {
  const classes = useStyles();

  const [inputValue, setInputValue] = useState('');
  const [active, setActive] = useState(false);
  const [disabled, setDisabled] = useState(false);

  const [dialogOpen, setDialogOpen] = useState(false);

  const dirty = !!filterValue.length;

  const addSmiles = smiles => {
    const newFilterValue = new Set([...filterValue, ...smiles]);
    setFilter([...newFilterValue]);
  };

  const { mutate: canonicalize } = useCanonicalizeSmiles(
    () => {
      setDisabled(true);
    },
    smiles => {
      if (!!smiles) {
        addSmiles(smiles);
        setInputValue('');
      }
      setDisabled(false);
    }
  );

  const addSmilesFromInput = () => {
    if (!!inputValue) {
      const smiles = inputValue.split(';').map(val => val.trim());
      const formData = new FormData();
      smiles.forEach(smile => formData.append('smiles', smile));
      canonicalize({ data: formData });
    }
  };

  return (
    <>
      <div className={classes.root}>
        <TextField
          id={id}
          label={label}
          variant="outlined"
          fullWidth
          value={inputValue}
          disabled={disabled}
          onChange={event => setInputValue(event.target.value)}
          onFocus={() => setActive(true)}
          onBlur={() => setActive(false)}
          InputProps={{
            className: classes.input,
            classes: { input: classes.inputField },
            startAdornment:
              !!filterValue.length &&
              filterValue.map(value => (
                <Chip
                  className={classes.chip}
                  key={value}
                  label={value}
                  tabIndex={-1}
                  disabled={disabled}
                  onDelete={() => {
                    const newFilterValue = [...filterValue];
                    const index = filterValue.findIndex(val => value === val);
                    newFilterValue.splice(index, 1);
                    setFilter(newFilterValue);
                  }}
                />
              )),
            endAdornment: (
              <div className={classes.endAdornment}>
                <IconButton
                  className={classNames(
                    classes.button,
                    classes.clear,
                    active && dirty && !disabled && classes.clearActive,
                    dirty && !disabled && 'smiles-filter-dirty'
                  )}
                  size="small"
                  title="Clear"
                  aria-label="Clear"
                  disabled={disabled}
                  onClick={() => setFilter([])}
                >
                  <Clear fontSize="small" />
                </IconButton>
                <IconButton
                  className={classes.button}
                  size="small"
                  title="Add smiles"
                  aria-label="Add smiles"
                  disabled={disabled}
                  onClick={addSmilesFromInput}
                >
                  <AddCircle fontSize="small" />
                </IconButton>
              </div>
            )
          }}
          inputProps={{
            autoComplete: 'off'
          }}
          onKeyPress={event => {
            if (event.key === 'Enter') {
              addSmilesFromInput();
            }
          }}
        />
        <Tooltip title="Upload smiles from file">
          <IconButton disabled={disabled} onClick={() => setDialogOpen(true)}>
            <CloudUpload />
          </IconButton>
        </Tooltip>
      </div>

      <CanonicalizeSmilesDialog
        open={dialogOpen}
        onClose={() => setDialogOpen(false)}
        onCanonicalizeStart={() => setDisabled(true)}
        onCanonicalizeEnd={smiles => {
          if (smiles) {
            addSmiles(smiles);
          }
          setDisabled(false);
        }}
      />
    </>
  );
};
