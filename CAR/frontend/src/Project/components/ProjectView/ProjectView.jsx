import { makeStyles } from '@material-ui/core';
import { MethodStepAccordion } from '../MethodStepAccordion';
import { useGetCategorizedMethodsBySteps } from './hooks/useGetCategorizedMethodsBySteps';

const useStyles = makeStyles((theme) => ({
  root: {
    padding: theme.spacing(2),
  },
}));

export const ProjectView = ({ projectId }) => {
  const classes = useStyles();

  const sortedMethods = useGetCategorizedMethodsBySteps(projectId);

  if (!projectId) {
    return null;
  }

  return (
    <main className={classes.root}>
      {Object.entries(sortedMethods).map(([noSteps, methods], index) => {
        return (
          <MethodStepAccordion
            key={noSteps}
            noSteps={Number(noSteps)}
            open={index === 0}
            methods={methods}
          />
        );
      })}
    </main>
  );
};
