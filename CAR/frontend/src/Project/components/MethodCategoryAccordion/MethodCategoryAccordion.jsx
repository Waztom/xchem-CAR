import {
  Accordion,
  AccordionSummary,
  colors,
  makeStyles,
} from '@material-ui/core';
import { ExpandMore } from '@material-ui/icons';
import { ImSad, ImSmile } from 'react-icons/im';
import { IoFootsteps } from 'react-icons/io5';

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
        className={classes.summary}
        classes={{
          content: classes.content,
        }}
        expandIcon={<ExpandMore />}
      >
        {successArray.map((success, index) => {
          return (
            <>
              <IoFootsteps className={classes.icon} />
              {success ? (
                <ImSmile key={index} className={classes.icon} />
              ) : (
                <ImSad key={index} className={classes.icon} />
              )}
            </>
          );
        })}
        <CategoryIcon className={classes.icon} />
      </AccordionSummary>
    </Accordion>
  );
};
