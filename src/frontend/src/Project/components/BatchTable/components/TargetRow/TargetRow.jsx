import React from 'react';
import { colors, TableCell, TableRow, Typography } from '@material-ui/core';
import { makeStyles } from '@material-ui/styles';
import { ExpandMore } from '@material-ui/icons';
import { useBatchContext } from '../../../../hooks/useBatchContext';
import { setRowsExpanded } from '../../../../../common/stores/batchesTableStateStore';

const useStyles = makeStyles(theme => ({
  root: {
    gridTemplateColumns: '40px repeat(2, auto) 1fr',
    justifyContent: 'flex-start',
    backgroundColor: colors.grey[100],
    cursor: 'pointer'
  },
  name: {
    width: '100%',
    fontWeight: 500
  },
  image: {
    mixBlendMode: 'multiply'
  },
  icon: {
    color: theme.palette.action.active,
    justifySelf: 'flex-end',
    transform: ({ expanded }) => `rotate(${expanded ? 180 : 0}deg)`
  }
}));

export const TargetRow = ({ row }) => {
  const classes = useStyles({ expanded: row.isExpanded });

  const batch = useBatchContext();

  const target = row.original;

  const selectionCell = row.cells.find(cell => cell.column.id === 'selection');

  return (
    <>
      <TableRow
        className={classes.root}
        onClick={() => {
          row.toggleRowExpanded();
          setRowsExpanded(batch.id, [row], !row.isExpanded);
        }}
      >
        <TableCell {...selectionCell.getCellProps()}>{selectionCell.render('Cell')}</TableCell>

        <TableCell>
          <img className={classes.image} src={target.image} width={120} height={60} alt={target.name} />
        </TableCell>

        <TableCell>
          <Typography className={classes.name} component="h3" noWrap>
            {target.name}
          </Typography>
        </TableCell>

        <TableCell>
          <ExpandMore className={classes.icon} />
        </TableCell>
      </TableRow>
    </>
  );
};
