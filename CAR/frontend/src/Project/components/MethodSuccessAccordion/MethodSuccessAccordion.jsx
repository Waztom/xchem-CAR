import {
  Accordion,
  AccordionDetails,
  AccordionSummary,
  colors,
  Divider,
  List,
  ListItem,
  makeStyles,
  Typography,
} from '@material-ui/core';
import { Cancel, ExpandMore, FindInPage } from '@material-ui/icons';
import { ImSad, ImSmile } from 'react-icons/im';
import { IoFootsteps } from 'react-icons/io5';
import { FaRegEdit, FaFlask } from 'react-icons/fa';
import { Fragment, useEffect } from 'react';
import axios from 'axios';
import { MethodCategoryAccordion } from '../MethodCategoryAccordion/MethodCategoryAccordion';

const useStyles = makeStyles((theme) => ({
  root: {
    width: '100%',
    boxShadow: 'none',
  },
  summary: {
    backgroundColor: colors.grey[300], // Might be removed
    position: 'sticky',
    top: 64,
    zIndex: 2,
  },
  content: {
    display: 'grid',
    gridTemplateColumns: 'minmax(max-content, 45%) minmax(max-content, 1fr)',
    '& > div': {
      display: 'flex',
      gap: theme.spacing(),
    },
  },
  details: {
    padding: 0,
  },
  categoryInfo: {
    display: 'flex',
    alignItems: 'center',
    gap: theme.spacing(1 / 2),
  },
  icon: {
    width: 24,
    height: 24,
  },
  list: {
    width: '100%',
    '& > li': {
      padding: 0,
    },
  },
}));

const temporaryData = [
  { CategoryIcon: FindInPage, value: 300 },
  { CategoryIcon: FaRegEdit, value: 10 },
  { CategoryIcon: FaFlask, value: 400 },
  { CategoryIcon: Cancel, value: 10 },
];

export const MethodSuccessAccordion = ({ successArray }) => {
  const classes = useStyles();

  useEffect(() => {
    async function fetchData() {
      const response = await axios.get(`/api/reactions/`);
      console.log(response);
    }
    //fetchData();
  }, []);

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
        <div>
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
        </div>
        <div>
          {temporaryData.map(({ value, CategoryIcon }, index) => {
            return (
              <div key={index} className={classes.categoryInfo}>
                <Typography>{value}</Typography>
                <CategoryIcon className={classes.icon} />
              </div>
            );
          })}
        </div>
      </AccordionSummary>
      <AccordionDetails className={classes.details}>
        <List className={classes.list} disablePadding>
          {temporaryData.map(({ CategoryIcon }, index) => {
            return (
              <Fragment key={index}>
                <ListItem disableGutters>
                  <MethodCategoryAccordion
                    successArray={successArray}
                    CategoryIcon={CategoryIcon}
                  />
                </ListItem>
                {!!(index < temporaryData.length - 1) && <Divider />}
              </Fragment>
            );
          })}
        </List>
      </AccordionDetails>
    </Accordion>
  );
};
