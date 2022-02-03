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
import { IconComponent } from '../../../common/components/IconComponent';
import { MethodSuccessAccordion } from '../MethodSuccessAccordion';
import { useCategorizeMethodsDataBySuccessRate } from './hooks/useCategorizeMethodsDataBySuccessRate';
import { useGetMethodsReactions } from './hooks/useGetMethodsReactions';

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

  const { methodsData } = useGetMethodsReactions(methodsWithTarget);
  const categorizedMethodsData =
    useCategorizeMethodsDataBySuccessRate(methodsData);

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
          <IconComponent key={index} Component={IoFootsteps} />
        ))}
      </AccordionSummary>
      <AccordionDetails className={classes.details}>
        <List className={classes.list} disablePadding>
          {Object.entries(categorizedMethodsData)
            .sort(([keyA], [keyB]) => keyB.localeCompare(keyA))
            .map(([successString, methodData], index) => {
              return (
                <Fragment key={successString}>
                  <ListItem disableGutters>
                    <MethodSuccessAccordion
                      noSteps={noSteps}
                      successString={successString}
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
