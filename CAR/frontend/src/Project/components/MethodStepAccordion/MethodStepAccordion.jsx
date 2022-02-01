import {
  Accordion,
  AccordionDetails,
  AccordionSummary,
  Divider,
  List,
  ListItem,
  makeStyles,
} from '@material-ui/core';
import { ExpandMore } from '@material-ui/icons';
import { Fragment, useLayoutEffect, useState } from 'react';
import { IoFootsteps } from 'react-icons/io5';
import { MethodSuccessAccordion } from '../MethodSuccessAccordion';
import { useCategorizeMethodReactionsBySuccess } from './hooks/useCategorizeMethodReactionsBySuccess';
import { useGetMethodReactions } from './hooks/useGetReactions';

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
  icon: {
    width: 24,
    height: 24,
  },
  details: {
    padding: 0,
  },
  list: {
    width: '100%',
    '& > li': {
      padding: 0,
    },
  },
}));

export const MethodStepAccordion = ({ noSteps, open, methods }) => {
  const classes = useStyles();

  const [expanded, setExpanded] = useState(open);

  const methodReactions = useGetMethodReactions(methods);
  const categorizedMethodReactions =
    useCategorizeMethodReactionsBySuccess(methodReactions);

  useLayoutEffect(() => {
    setExpanded(expanded);
  }, [open]);

  return (
    <Accordion
      expanded={expanded}
      onChange={(_, expanded) => setExpanded(expanded)}
    >
      <AccordionSummary
        className={classes.summary}
        expandIcon={<ExpandMore className={classes.collapseIcon} />}
      >
        {new Array(noSteps).fill(0).map((_, index) => (
          <IoFootsteps key={index} className={classes.icon} />
        ))}
      </AccordionSummary>
      <AccordionDetails className={classes.details}>
        <List className={classes.list} disablePadding>
          {Object.entries(categorizedMethodReactions)
            .reverse()
            .map(([noSuccesses, methodReactions], index) => {
              return (
                <Fragment key={noSuccesses}>
                  <ListItem disableGutters>
                    <MethodSuccessAccordion
                      noSteps={noSteps}
                      noSuccesses={Number(noSuccesses)}
                      methodReactions={methodReactions}
                    />
                  </ListItem>
                  {!!(index < methodReactions.length - 1) && <Divider />}
                </Fragment>
              );
            })}
        </List>
      </AccordionDetails>
    </Accordion>
  );
};
