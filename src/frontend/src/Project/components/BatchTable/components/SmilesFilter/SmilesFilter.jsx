import React, { useState } from 'react';
import { Chip, IconButton, TextField, Tooltip } from '@mui/material';
import { styled } from '@mui/material/styles';
import { AddCircle, Clear, CloudUpload } from '@mui/icons-material';
import { SuspenseWithBoundary } from '../../../../../common/components/SuspenseWithBoundary';
import { CanonicalizeSmilesDialog } from '../../../CanonicalizeSmilesDialog';
import { useCanonicalizeSmiles } from '../../../../hooks/useCanonicalizeSmiles';

const FilterRoot = styled('div')(({ theme }) => ({
  display: 'flex',
  alignItems: 'center',
  gap: theme.spacing()
}));

const StyledTextField = styled(TextField)(({ theme }) => ({
  '& .MuiInputBase-root': {
    padding: 9,
    paddingRight: 65,
    flexWrap: 'wrap',
    '&:hover .smiles-filter-dirty': {
      visibility: 'visible'
    }
  },
  '& .MuiInputBase-input': {
    flexGrow: 1,
    minWidth: 30,
    width: 0,
    padding: '9.5px 4px'
  }
}));

const StyledChip = styled(Chip)(({ theme }) => ({
  margin: 3,
  maxWidth: 240
}));

const EndAdornment = styled('div')(({ theme }) => ({
  position: 'absolute',
  right: 9,
  top: 'calc(50% - 14px)',
  display: 'flex',
  gap: theme.spacing(1/2)
}));

const ActionButton = styled(IconButton, {
  shouldForwardProp: prop => prop !== 'isVisible'
})(({ theme, isVisible }) => ({
  padding: 2,
  marginRight: -2,
  visibility: isVisible ? 'visible' : 'hidden'
}));

const SmilesFilterContent = ({ id, label, filterValue = [], setFilter }) => {
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
    () => setDisabled(true),
    smiles => {
      if (smiles) {
        addSmiles(smiles);
        setInputValue('');
      }
      setDisabled(false);
    }
  );

  const addSmilesFromInput = () => {
    if (inputValue) {
      const smiles = inputValue.split(';').map(val => val.trim());
      const formData = new FormData();
      smiles.forEach(smile => formData.append('smiles', smile));
      canonicalize({ data: formData });
    }
  };

  return (
    <FilterRoot>
      <StyledTextField
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
          startAdornment: filterValue.map(value => (
            <StyledChip
              key={value}
              label={value}
              tabIndex={-1}
              disabled={disabled}
              onDelete={() => {
                setFilter(filterValue.filter(val => val !== value));
              }}
            />
          )),
          endAdornment: (
            <EndAdornment>
              <ActionButton
                size="small"
                title="Clear"
                isVisible={active && dirty && !disabled}
                disabled={disabled}
                onClick={() => setFilter([])}
              >
                <Clear fontSize="small" />
              </ActionButton>
              <ActionButton
                size="small"
                title="Add SMILES"
                isVisible={true}
                disabled={disabled}
                onClick={addSmilesFromInput}
              >
                <AddCircle fontSize="small" />
              </ActionButton>
            </EndAdornment>
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
      <Tooltip title="Upload SMILES from file">
        <IconButton 
          disabled={disabled} 
          onClick={() => setDialogOpen(true)}
          size="large"
        >
          <CloudUpload />
        </IconButton>
      </Tooltip>

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
    </FilterRoot>
  );
};

export const SmilesFilter = (props) => (
  <SuspenseWithBoundary>
    <SmilesFilterContent {...props} />
  </SuspenseWithBoundary>
);

SmilesFilter.displayName = 'SmilesFilter';
