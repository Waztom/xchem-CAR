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
import axios from 'axios';
import { Fragment, useEffect, useLayoutEffect, useState } from 'react';
import { IoFootsteps } from 'react-icons/io5';
import { MethodSuccessAccordion } from '../MethodSuccessAccordion';
import { colors } from '@material-ui/core';

const useStyles = makeStyles((theme) => ({
  summary: {
    backgroundColor: theme.palette.primary.main,
    color: theme.palette.white,
  },
  collapseIcon: {
    color: theme.palette.white,
  },
  icon: {
    width: 24,
    height: 24,
  },
  accordionDetails: {
    padding: 0,
  },
  list: {
    width: '100%',
    '& > li': {
      backgroundColor: colors.grey[300], // Might be removed
      padding: 0,
    },
  },
}));

/**
 * Creates a permutation of success rates of steps in a triangular pattern, e.g.:
 * [
 *   [true,   true,   true],
 *   [true,   true,   false],
 *   [true,   false,  false],
 *   [false,  false,  false],
 * ]
 * @param {number} noSteps Number of steps
 * @returns Permutation of true / false values.
 */
const getSuccessRatePermutation = (noSteps) => {
  const successPermutation = [];
  for (let i = 0; i <= noSteps; i++) {
    successPermutation[i] = [];
    for (let j = 0; j < noSteps; j++) {
      successPermutation[i][j] = j < noSteps - i;
    }
  }
  return successPermutation;
};

export const MethodAccordion = ({ noSteps, open }) => {
  const classes = useStyles();

  const [expanded, setExpanded] = useState(open);

  const successPermutation = getSuccessRatePermutation(noSteps);

  useLayoutEffect(() => {
    setExpanded(expanded);
  }, [open]);

  useEffect(() => {
    async function fetchData() {
      const response = await axios.get(`/api/methods/?nosteps=${noSteps}`);
      console.log(response);
    }
    fetchData();
  }, []);

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
      <AccordionDetails className={classes.accordionDetails}>
        <List className={classes.list} disablePadding>
          {successPermutation.map((successArray, index) => {
            return (
              <Fragment key={index}>
                <ListItem disableGutters>
                  <MethodSuccessAccordion successArray={successArray} />
                </ListItem>
                {!!(index < successPermutation.length - 1) && <Divider />}
              </Fragment>
            );
          })}
        </List>
      </AccordionDetails>
    </Accordion>
  );
};
