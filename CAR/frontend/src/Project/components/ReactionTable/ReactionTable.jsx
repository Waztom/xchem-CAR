import {
  Checkbox,
  IconButton,
  makeStyles,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  TableSortLabel,
  Typography,
} from '@material-ui/core';
import { useMemo } from 'react';
import { useTable, useSortBy } from 'react-table';
import { SiMoleculer } from 'react-icons/si';
import { IconComponent } from '../../../common/components/IconComponent';
import { FaFlask } from 'react-icons/fa';
import { GiMoneyStack } from 'react-icons/gi';
import { IoFootsteps } from 'react-icons/io5';
import { ImSad, ImSmile } from 'react-icons/im';

const useStyles = makeStyles((theme) => ({
  table: {
    width: 'auto',
  },
  cell: {
    textAlign: 'center',
  },
  sortCell: {
    display: 'flex',
    justifyContent: 'center',
    alignItems: 'center',
    gap: theme.spacing(1 / 2),
  },
}));

export const ReactionTable = ({ noSteps, methodData }) => {
  const classes = useStyles();

  console.log(methodData);

  const columns = useMemo(() => {
    return [
      {
        accessor: 'target.image',
        disableSortBy: true,
        Header: () => {
          return <IconComponent Component={SiMoleculer} />;
        },
        Cell: ({ value }) => {
          return <img src={value} height={60} />;
        },
      },
      {
        accessor: 'method.estimatecost',
        Header: () => {
          return <IconComponent Component={GiMoneyStack} />;
        },
      },
      {
        accessor: 'method.synthesise',
        disableSortBy: true,
        Header: () => {
          return <IconComponent Component={FaFlask} />;
        },
        Cell: ({ value }) => {
          return <Checkbox checked={value} />;
        },
      },
      ...new Array(noSteps).fill(0).map((_, index) => {
        return {
          accessor: `reactions[${index}].reactionimage`,
          disableSortBy: true,
          Header: () => {
            return (
              <div className={classes.sortCell}>
                <IconComponent Component={IoFootsteps} />
                <Typography component="span">{index + 1}</Typography>
              </div>
            );
          },
          Cell: ({ value }) => {
            return (
              <>
                <img src={value} height={60} />
                <IconButton size="small">
                  <IconComponent Component={ImSmile} />
                </IconButton>
                <IconButton size="small">
                  <IconComponent Component={ImSad} />
                </IconButton>
              </>
            );
          },
        };
      }),
    ];
  }, [noSteps]);

  const { getTableProps, getTableBodyProps, headerGroups, prepareRow, rows } =
    useTable(
      {
        columns,
        data: methodData,
      },
      useSortBy
    );

  return (
    <Table className={classes.table} {...getTableProps()}>
      <TableHead>
        {headerGroups.map((headerGroup) => (
          <TableRow {...headerGroup.getHeaderGroupProps()}>
            {headerGroup.headers.map((column) => (
              <TableCell
                className={classes.cell}
                {...column.getHeaderProps(
                  column.canSort ? column.getSortByToggleProps() : undefined
                )}
              >
                {column.render('Header')}
                {column.canSort && (
                  <TableSortLabel
                    active={column.isSorted}
                    // react-table has a unsorted state which is not treated here
                    direction={column.isSortedDesc ? 'desc' : 'asc'}
                  />
                )}
              </TableCell>
            ))}
          </TableRow>
        ))}
      </TableHead>
      <TableBody {...getTableBodyProps()}>
        {rows.map((row, i) => {
          prepareRow(row);
          return (
            <TableRow {...row.getRowProps()}>
              {row.cells.map((cell) => (
                <TableCell className={classes.cell} {...cell.getCellProps()}>
                  {cell.render('Cell')}
                </TableCell>
              ))}
            </TableRow>
          );
        })}
      </TableBody>
    </Table>
  );
};
