import { makeStyles } from '@material-ui/core';
import { MethodAccordion } from '../MethodAccordion';

const useStyles = makeStyles((theme) => ({
  root: {
    padding: theme.spacing(2),
  },
}));

const methodAccordions = [
  { noSteps: 1, open: true },
  { noSteps: 2, open: false },
  { noSteps: 3, open: false },
];

export const ProjectView = ({ projectId }) => {
  const classes = useStyles();

  if (!projectId) {
    return null;
  }

  return (
    <main className={classes.root}>
      {methodAccordions.map(({ noSteps, open }) => {
        return <MethodAccordion key={noSteps} noSteps={noSteps} open={open} />;
      })}
    </main>
  );
};
