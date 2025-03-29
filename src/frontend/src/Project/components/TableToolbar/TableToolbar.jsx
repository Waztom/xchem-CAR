import React, { Fragment, useCallback, useState } from 'react';
import { Accordion, AccordionDetails, AccordionSummary, makeStyles, Typography } from '@material-ui/core';
import { useBatchesTableStateStore } from '../../../common/stores/batchesTableStateStore';
import { useBatchContext } from '../../hooks/useBatchContext';
import { ToolbarSection } from './components/ToolbarSection';
import { ExpandMore } from '@material-ui/icons';
import { ButtonSection } from './components/ButtonSection';

const useStyles = makeStyles(theme => ({
  root: {
    margin: '0 !important'
  },
  heading: {
    fontWeight: 500
  },
  details: {
    display: 'flex',
    gap: theme.spacing(2),
    overflow: 'auto',
    padding: `${theme.spacing()}px ${theme.spacing(2)}px`,
    // Filters column
    '& > :nth-child(2)': {
      flex: '1 0 320px'
    }
  },
  firstColumn: {
    flex: '0 0 260px'
  },
  filterColumn: {
    display: 'grid',
    gap: theme.spacing(),
    width: '100%'
  }
}));

export const TableToolbar = ({ tableInstance }) => {
  const classes = useStyles();

  const batch = useBatchContext();

  const { columns } = tableInstance;

  const selectedMethodsIds = useBatchesTableStateStore(
    useCallback(
      state =>
        Object.entries(state.selected[batch.id] || {})
          .filter(([_, value]) => value)
          // Method row ids are in a form targetId.methodId
          .map(([key]) => Number(key.split('.')[1])),
      [batch.id]
    )
  );

  const [accordionOpen, setAccordionOpen] = useState(false);

  return (
    <Accordion
      expanded={accordionOpen}
      className={classes.root}
      elevation={0}
      onChange={(event, expanded) => setAccordionOpen(expanded)}
    >
      <AccordionSummary expandIcon={<ExpandMore />}>
        <Typography className={classes.heading}>
          {accordionOpen ? 'Hide' : 'Show'} table summary, actions and filters
        </Typography>
      </AccordionSummary>
      <AccordionDetails className={classes.details}>
        <div className={classes.firstColumn}>
          <ToolbarSection title="Summary">
            <Typography>Selected methods: {selectedMethodsIds.length}</Typography>
          </ToolbarSection>
          <ButtonSection tableInstance={tableInstance} selectedMethodsIds={selectedMethodsIds} />
        </div>
        <ToolbarSection title="Filters">
          {columns
            .filter(column => column.canFilter)
            .sort((a, b) => a.filterOrder - b.filterOrder)
            .map(column => (
              <Fragment key={column.id}>{column.render('Filter')}</Fragment>
            ))}
        </ToolbarSection>
      </AccordionDetails>
    </Accordion>
  );
};
