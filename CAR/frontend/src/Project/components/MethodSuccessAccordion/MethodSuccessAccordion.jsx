import {
  Accordion,
  AccordionSummary,
  List,
  ListItem,
  makeStyles,
  Typography,
} from '@material-ui/core';
import { Cancel, ExpandMore, FindInPage } from '@material-ui/icons';
import { ImSad, ImSmile } from 'react-icons/im';
import { IoFootsteps } from 'react-icons/io5';
import { FaRegEdit, FaFlask } from 'react-icons/fa';

const useStyles = makeStyles((theme) => ({
  root: {
    width: '100%',
    boxShadow: 'none',
  },
  accordionSummary: {
    display: 'grid',
    gridTemplateColumns: 'minmax(max-content, 45%) minmax(max-content, 1fr)',
    '& > div': {
      display: 'flex',
      gap: theme.spacing(),
      alignItems: 'center',
    },
  },
  icon: {
    width: 24,
    height: 24,
  },
}));

export const MethodSuccessAccordion = ({ successArray }) => {
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
        aria-controls={`method${0}-succes-content`}
        id={`method${0}-succes-header`}
      >
        <div>
          <IoFootsteps className={classes.icon} />
          {successArray.map((success, index) => {
            return success ? (
              <ImSmile key={index} className={classes.icon} />
            ) : (
              <ImSad key={index} className={classes.icon} />
            );
          })}
        </div>
        <div>
          <Typography>300</Typography>
          <FindInPage className={classes.icon} />
          <Typography>10</Typography>
          <FaRegEdit className={classes.icon} />
          <Typography>400</Typography>
          <FaFlask className={classes.icon} />
          <Typography>10</Typography>
          <Cancel className={classes.icon} />
        </div>
      </AccordionSummary>
    </Accordion>
  );
};
