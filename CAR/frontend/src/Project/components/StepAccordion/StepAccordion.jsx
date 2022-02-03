import {
  Accordion,
  AccordionDetails,
  AccordionSummary,
  makeStyles,
} from '@material-ui/core';
import { ExpandMore } from '@material-ui/icons';
import { useLayoutEffect, useState } from 'react';
import { IoFootsteps } from 'react-icons/io5';
import { IconComponent } from '../../../common/components/IconComponent';
import { SuccessRateAccordionList } from './components/SuccessRateAccordionList';

const useStyles = makeStyles((theme) => ({
  summary: {
    backgroundColor: theme.palette.primary.main,
    color: theme.palette.white,
    position: 'sticky',
    top: 0,
    zIndex: 3,
  },
  collapseIcon: {
    color: theme.palette.white,
  },
  details: {
    padding: 0,
  },
}));

export const StepAccordion = ({ noSteps, open, methodsWithTarget }) => {
  const classes = useStyles();

  const [expanded, setExpanded] = useState(open);

  useLayoutEffect(() => {
    setExpanded(expanded);
  }, [open]);

  return (
    <Accordion
      expanded={expanded}
      onChange={(_, expanded) => setExpanded(expanded)}
      TransitionProps={{ unmountOnExit: true }} // Performance
    >
      <AccordionSummary
        className={classes.summary}
        expandIcon={<ExpandMore className={classes.collapseIcon} />}
      >
        {new Array(noSteps).fill(0).map((_, index) => (
          <IconComponent key={index} Component={IoFootsteps} />
        ))}
      </AccordionSummary>
      <AccordionDetails className={classes.details}>
        <SuccessRateAccordionList
          noSteps={noSteps}
          methodsWithTarget={methodsWithTarget}
        />
      </AccordionDetails>
    </Accordion>
  );
};
