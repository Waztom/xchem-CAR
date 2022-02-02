import {
  Accordion,
  AccordionDetails,
  AccordionSummary,
  colors,
  makeStyles,
} from '@material-ui/core';
import { ExpandMore } from '@material-ui/icons';
import { Fragment } from 'react';
import { ImSad, ImSmile } from 'react-icons/im';
import { IoFootsteps } from 'react-icons/io5';
import { ReactionTable } from '../ReactionTable';

const useStyles = makeStyles((theme) => ({
  root: {
    width: '100%',
    boxShadow: 'none',
  },
  summary: {
    backgroundColor: colors.grey[100], // Might be removed
    position: 'sticky',
    top: 128,
    zIndex: 1,
  },
  content: {
    display: 'flex',
    gap: theme.spacing(),
  },
  details: {
    padding: 0,
  },
  icon: {
    width: 24,
    height: 24,
  },
}));

export const MethodCategoryAccordion = ({
  noSteps,
  noSuccesses,
  methodData,
  CategoryIcon,
}) => {
  const classes = useStyles();

  return (
    <Accordion
      className={classes.root}
      TransitionProps={{ unmountOnExit: true }} // Performance
    >
      <AccordionSummary
        className={classes.summary}
        classes={{
          content: classes.content,
        }}
        expandIcon={<ExpandMore />}
      >
        {new Array(noSuccesses).fill(0).map((_, index) => {
          return (
            <Fragment key={index}>
              <IoFootsteps className={classes.icon} />
              <ImSmile key={index} className={classes.icon} />
            </Fragment>
          );
        })}
        {new Array(noSteps - noSuccesses).fill(0).map((_, index) => {
          return (
            <Fragment key={index}>
              <IoFootsteps className={classes.icon} />
              <ImSad key={index} className={classes.icon} />
            </Fragment>
          );
        })}
        <CategoryIcon className={classes.icon} />
      </AccordionSummary>
      <AccordionDetails className={classes.details}>
        <ReactionTable noSteps={noSteps} methodData={methodData} />
      </AccordionDetails>
    </Accordion>
  );
};
