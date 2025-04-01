import React from 'react';
import {
  colors,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  TableSortLabel,
  Tooltip,
  Typography,
  Box
} from '@mui/material';
import { styled } from '@mui/material/styles';
import { DialogSection } from '../../../../../common/components/DialogSection';
import { DialogSectionHeading } from '../../../../../common/components/DialogSectionHeading';
import { useReactantProductTableColumns } from './hooks/useReactantProductTableColumns';
import { useTable, useSortBy } from 'react-table';
import { SuspenseWithBoundary } from '../../../../../common/components/SuspenseWithBoundary';
import PubChem from '../../../../../assets/pubchem.svg';
import PubChemSafety from '../../../../../assets/pubchem-safety.svg';

const StyledTable = styled(Table)(({ theme }) => ({
  '& tr': {
    display: 'grid',
    alignItems: 'stretch',
    gap: `0 ${theme.spacing()}px`,
    borderBottom: `1px solid ${theme.palette.divider}`,
    paddingLeft: theme.spacing(2),
    paddingRight: theme.spacing(2),
    gridTemplateColumns: '160px 1fr 105px repeat(2, 95px)'
  },
  '& th, td': {
    display: 'grid',
    placeItems: 'center start',
    border: 0,
    padding: 0
  }
}));

const FlexCell = styled('div')(({ theme }) => ({
  display: 'flex',
  alignItems: 'center',
  gap: theme.spacing(1/2)
}));

const StyledSortLabel = styled(TableSortLabel)(({ theme, active }) => ({
  '& svg': {
    color: !active ? `${colors.grey[400]} !important` : undefined
  }
}));

const SectionsGrid = styled('div')(({ theme }) => ({
  display: 'grid',
  gridTemplateColumns: 'repeat(3, 1fr)'
}));

const ReactantSectionContent = ({ reactant, index }) => {
  const reactantpubcheminfo = reactant.reactantpubcheminfo || {};
  const columns = useReactantProductTableColumns();
  const tableInstance = useTable(
    {
      columns,
      data: reactant.catalogentries
    },
    useSortBy
  );

  const { getTableProps, headerGroups, getTableBodyProps, rows, prepareRow } = tableInstance;

  return (
    <DialogSection>
      <DialogSectionHeading>Reactant {index + 1}</DialogSectionHeading>
      <SectionsGrid>
        <Box>
          <Typography>
            Smiles: <strong>{reactant.smiles}</strong>
          </Typography>
          {reactantpubcheminfo.cas && (
            <Typography>
              CAS number: <strong>{reactantpubcheminfo.cas}</strong>
            </Typography>
          )}
        </Box>

        {!!reactantpubcheminfo.summaryurl && (
          <Box>
            <Tooltip title="PubChem compound summary">
              <a href={reactantpubcheminfo.summaryurl} target="_blank" rel="noreferrer">
                <img src={PubChem} height={50} alt="PubChem" />
              </a>
            </Tooltip>
          </Box>
        )}
        {!!reactantpubcheminfo.lcssurl && (
          <Box>
            <Tooltip title="Laboratory chemical safety summary">
              <a href={reactantpubcheminfo.lcssurl} target="_blank" rel="noreferrer">
                <img src={PubChemSafety} height={50} alt="PubChem Safety" />
              </a>
            </Tooltip>
          </Box>
        )}
      </SectionsGrid>

      {!!reactant.catalogentries.length && (
        <>
          <Typography>Catalog information: </Typography>
          <StyledTable {...getTableProps()}>
            <TableHead>
              {headerGroups.map(headerGroup => (
                <TableRow {...headerGroup.getHeaderGroupProps()}>
                  {headerGroup.headers.map(column => {
                    if (column.canSort) {
                      const { title, ...rest } = column.getSortByToggleProps();
                      return (
                        <Tooltip title={`Sort by ${column.sortLabel}`} {...column.getHeaderProps()}>
                          <TableCell {...rest}>
                            <FlexCell>
                              {column.render('Header')}
                              <StyledSortLabel
                                active={true}
                                direction={column.isSortedDesc ? 'desc' : 'asc'}
                              />
                            </FlexCell>
                          </TableCell>
                        </Tooltip>
                      );
                    }
                    return <TableCell {...column.getHeaderProps()}>{column.render('Header')}</TableCell>;
                  })}
                </TableRow>
              ))}
            </TableHead>
            <TableBody {...getTableBodyProps()}>
              {rows.map(row => {
                prepareRow(row);
                return (
                  <TableRow {...row.getRowProps()}>
                    {row.cells.map(cell => (
                      <TableCell {...cell.getCellProps()}>{cell.render('Cell')}</TableCell>
                    ))}
                  </TableRow>
                );
              })}
            </TableBody>
          </StyledTable>
        </>
      )}
    </DialogSection>
  );
};

export const ReactantSection = (props) => (
  <SuspenseWithBoundary>
    <ReactantSectionContent {...props} />
  </SuspenseWithBoundary>
);

ReactantSection.displayName = 'ReactantSection';
