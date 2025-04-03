import React, { Fragment, useCallback, useState } from 'react';
import { Accordion, AccordionDetails, AccordionSummary, Typography } from '@mui/material';
import { styled } from '@mui/material/styles';
import { ExpandMore } from '@mui/icons-material';
import { SuspenseWithBoundary } from '../../../common/components/SuspenseWithBoundary';
import { useBatchesTableStateStore } from '../../../common/stores/batchesTableStateStore';
import { useBatchContext } from '../../hooks/useBatchContext';
import { ToolbarSection } from './components/ToolbarSection';
import { ButtonSection } from './components/ButtonSection';

const StyledAccordion = styled(Accordion)(({ theme }) => ({
  margin: '0 !important'
}));

const StyledAccordionDetails = styled(AccordionDetails)(({ theme }) => ({
  display: 'flex',
  gap: theme.spacing(2),
  overflow: 'auto',
  padding: `${theme.spacing()}px ${theme.spacing(2)}`,
  '& > *:nth-of-type(2)': {
    flex: '1 0 320px'
  }
}));

const FirstColumn = styled('div')(({ theme }) => ({
  flex: '0 0 260px'
}));

const FilterColumn = styled('div')(({ theme }) => ({
  display: 'grid',
  gap: theme.spacing(),
  width: '100%'
}));

const StyledHeading = styled(Typography)(({ theme }) => ({
  fontWeight: 500
}));

const TableToolbarContent = ({ tableInstance }) => {
  const batch = useBatchContext();
  const { columns } = tableInstance;
  const [accordionOpen, setAccordionOpen] = useState(false);

  const selectedMethodsIds = useBatchesTableStateStore(
    useCallback(
      state =>
        Object.entries(state.selected[batch.id] || {})
          .filter(([_, value]) => value)
          .map(([key]) => Number(key.split('.')[1])),
      [batch.id]
    )
  );

  return (
    <StyledAccordion
      expanded={accordionOpen}
      elevation={0}
      onChange={(_, expanded) => setAccordionOpen(expanded)}
    >
      <AccordionSummary expandIcon={<ExpandMore />}>
        <StyledHeading>
          {accordionOpen ? 'Hide' : 'Show'} table summary, actions and filters
        </StyledHeading>
      </AccordionSummary>
      <StyledAccordionDetails>
        <FirstColumn>
          <ToolbarSection title="Summary">
            <Typography>Selected methods: {selectedMethodsIds.length}</Typography>
          </ToolbarSection>
          <SuspenseWithBoundary>
            <ButtonSection 
              tableInstance={tableInstance} 
              selectedMethodsIds={selectedMethodsIds} 
            />
          </SuspenseWithBoundary>
        </FirstColumn>
        <ToolbarSection title="Filters">
          <FilterColumn>
            {columns
              .filter(column => column.canFilter)
              .sort((a, b) => a.filterOrder - b.filterOrder)
              .map(column => (
                <SuspenseWithBoundary key={column.id}>
                  {column.render('Filter')}
                </SuspenseWithBoundary>
              ))}
          </FilterColumn>
        </ToolbarSection>
      </StyledAccordionDetails>
    </StyledAccordion>
  );
};

export const TableToolbar = (props) => (
  <SuspenseWithBoundary>
    <TableToolbarContent {...props} />
  </SuspenseWithBoundary>
);

TableToolbar.displayName = 'TableToolbar';
