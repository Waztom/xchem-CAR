import React from 'react';
import { colors, TableCell, TableRow, Typography } from '@mui/material';
import { styled } from '@mui/material/styles';
import { ExpandMore } from '@mui/icons-material';
import { useBatchContext } from '../../../../hooks/useBatchContext';
import { setRowsExpanded } from '../../../../../common/stores/batchesTableStateStore';
import { SuspenseWithBoundary } from '../../../../../common/components/SuspenseWithBoundary';

const StyledTableRow = styled(TableRow)(({ theme }) => ({
  gridTemplateColumns: '40px repeat(2, auto) 1fr',
  justifyContent: 'flex-start',
  backgroundColor: colors.grey[100],
  cursor: 'pointer'
}));

const StyledName = styled(Typography)(() => ({
  width: '100%',
  fontWeight: 500
}));

const StyledImage = styled('img')(() => ({
  mixBlendMode: 'multiply'
}));

const StyledIcon = styled(ExpandMore, {
  shouldForwardProp: prop => prop !== 'expanded'
})(({ theme, expanded }) => ({
  color: theme.palette.action.active,
  justifySelf: 'flex-end',
  transform: `rotate(${expanded ? 180 : 0}deg)`,
  transition: theme.transitions.create('transform')
}));

const TargetRowContent = ({ row }) => {
  const batch = useBatchContext();
  const target = row.original;
  const selectionCell = row.cells.find(cell => cell.column.id === 'selection');
  const { key: selectionKey, ...selectionProps } = selectionCell.getCellProps();

  return (
    <StyledTableRow
      onClick={() => {
        row.toggleRowExpanded();
        setRowsExpanded(batch.id, [row], !row.isExpanded);
      }}
    >
      <TableCell key={selectionKey} {...selectionProps}>
        {selectionCell.render('Cell')}
      </TableCell>

      <TableCell>
        <StyledImage 
          src={target.image} 
          width={120} 
          height={60} 
          alt={target.name} 
        />
      </TableCell>

      <TableCell>
        <StyledName component="h3" noWrap>
          {target.name}
        </StyledName>
      </TableCell>

      <TableCell>
        <StyledIcon expanded={row.isExpanded} />
      </TableCell>
    </StyledTableRow>
  );
};

export const TargetRow = (props) => (
  <SuspenseWithBoundary>
    <TargetRowContent {...props} />
  </SuspenseWithBoundary>
);

TargetRow.displayName = 'TargetRow';
