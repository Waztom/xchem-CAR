import React from 'react';
import { Checkbox, Tooltip, Typography } from '@mui/material';
import { styled } from '@mui/material/styles';
import { TreeItem } from '@mui/x-tree-view/TreeItem';

const StyledTreeItem = styled(TreeItem)(({ theme }) => ({
  '& .MuiTreeItem-label': {
    display: 'flex',
    minWidth: 0,
    alignItems: 'center'
  },
  '& .MuiTreeItem-content': {
    '&.leaf': {
      cursor: 'default'
    }
  }
}));

const StyledTypography = styled(Typography)({
  flexGrow: 1
});

const StyledCheckbox = styled(Checkbox)({
  padding: 0
});

export const BatchSelectorItem = ({ batch, children, selected, onSelect }) => {
  return (
    <StyledTreeItem
      nodeId={String(batch.id)}
      className={!children.length ? 'leaf' : ''}
      label={
        <>
          <StyledTypography noWrap>
            {batch.batchtag}
          </StyledTypography>
          <Tooltip title={selected ? 'Deselect batch' : 'Select batch'}>
            <StyledCheckbox
              checked={selected}
              onClick={e => e.stopPropagation()}
              onChange={(_, checked) => onSelect(batch.id, checked)}
              inputProps={{ 'aria-label': batch.batchtag }}
            />
          </Tooltip>
        </>
      }
    >
      {children}
    </StyledTreeItem>
  );
};

BatchSelectorItem.displayName = 'BatchSelectorItem';
