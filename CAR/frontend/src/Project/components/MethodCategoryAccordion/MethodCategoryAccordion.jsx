import { Accordion, AccordionSummary, makeStyles } from '@material-ui/core';
import { ExpandMore } from '@material-ui/icons';
import { ImSad, ImSmile } from 'react-icons/im';
import { IoFootsteps } from 'react-icons/io5';

const useStyles = makeStyles((theme) => ({
  root: {
    width: '100%',
    boxShadow: 'none',
    backgroundColor: 'transparent', // Might be removed
  },
  accordionSummary: {
    display: 'flex',
    gap: theme.spacing(),
  },
  icon: {
    width: 24,
    height: 24,
  },
}));

export const MethodCategoryAccordion = ({ successArray, CategoryIcon }) => {
  const classes = useStyles();

  return (
    <Accordion
      className={classes.root}
      TransitionProps={{ unmountOnExit: true }} // Performance
    >
      <AccordionSummary
        classes={{
          content: classes.accordionSummary,
        }}
        expandIcon={<ExpandMore />}
      >
        <IoFootsteps className={classes.icon} />
        {successArray.map((success, index) => {
          return success ? (
            <ImSmile key={index} className={classes.icon} />
          ) : (
            <ImSad key={index} className={classes.icon} />
          );
        })}
        <CategoryIcon className={classes.icon} />
      </AccordionSummary>
    </Accordion>
  );
};
