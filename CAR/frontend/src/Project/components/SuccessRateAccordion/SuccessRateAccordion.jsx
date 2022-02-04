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
import { Fragment } from 'react';
import { CategoryAccordion } from '../CategoryAccordion';
import { IconComponent } from '../../../common/components/IconComponent';

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
  list: {
    width: '100%',
    '& > li': {
      padding: 0,
    },
  },
}));

const temporaryData = [
  { type: 'review', CategoryIcon: FindInPage, value: 300 },
  { type: 'edit', CategoryIcon: FaRegEdit, value: 10 },
  { type: 'synthesise', CategoryIcon: FaFlask, value: 400 },
  { type: 'ignore', CategoryIcon: Cancel, value: 10 },
];

export const SuccessRateAccordion = ({
  noSteps,
  successString,
  methodData,
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
        aria-controls={`successrate-accordion-${noSteps}-${successString}-content`}
        id={`successrate-accordion-${noSteps}-${successString}-header`}
      >
        <div>
          {[...successString].map((successChar, index) => {
            return (
              <Fragment key={index}>
                <IconComponent Component={IoFootsteps} />
                <IconComponent
                  Component={successChar === '1' ? ImSmile : ImSad}
                />
              </Fragment>
            );
          })}
        </div>
        <div>
          {temporaryData.map(({ value, CategoryIcon }, index) => {
            return (
              <div key={index} className={classes.categoryInfo}>
                <Typography>{value}</Typography>
                <IconComponent Component={CategoryIcon} />
              </div>
            );
          })}
        </div>
      </AccordionSummary>
      <AccordionDetails className={classes.details}>
        <List className={classes.list} disablePadding>
          {temporaryData.map((category, index) => {
            return (
              <Fragment key={index}>
                <ListItem disableGutters>
                  <CategoryAccordion
                    noSteps={noSteps}
                    successString={successString}
                    methodData={methodData}
                    category={category}
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
