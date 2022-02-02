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
import { useCategorizeMethodDataBySuccess } from './hooks/useCategorizeMethodDataBySuccess';
import { useGetMethodReactions } from './hooks/useGetMethodReactions';

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

export const MethodStepAccordion = ({ noSteps, open, methodsWithTarget }) => {
  const classes = useStyles();

  const [expanded, setExpanded] = useState(open);

  const methodsData = useGetMethodReactions(methodsWithTarget);
  const categorizedMethodsData =
    useCategorizeMethodDataBySuccess(methodsData);

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
          {Object.entries(categorizedMethodsData)
            .reverse()
            .map(([noSuccesses, methodData], index) => {
              return (
                <Fragment key={noSuccesses}>
                  <ListItem disableGutters>
                    <MethodSuccessAccordion
                      noSteps={noSteps}
                      noSuccesses={Number(noSuccesses)}
                      methodData={methodData}
                    />
                  </ListItem>
                  {!!(index < methodData.length - 1) && <Divider />}
                </Fragment>
              );
            })}
        </List>
      </AccordionDetails>
    </Accordion>
  );
};
