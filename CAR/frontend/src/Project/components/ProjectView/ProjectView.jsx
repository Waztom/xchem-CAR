import { makeStyles } from '@material-ui/core';
import { MethodStepAccordion } from '../MethodStepAccordion';
import { useCategorizeMethodsByNoSteps } from './hooks/useCategorizeMethodsByNoSteps';
import { useGetMethodsForTargets } from './hooks/useGetMethodsForTargets';
import { useGetTargets } from './hooks/useGetTargets';

const useStyles = makeStyles((theme) => ({
  root: {
    padding: theme.spacing(2),
  },
}));

export const ProjectView = () => {
  const classes = useStyles();

  const { data: targets } = useGetTargets();
  const { methodsWithTarget } = useGetMethodsForTargets(targets);
  const categorizedMethodsWithTarget =
    useCategorizeMethodsByNoSteps(methodsWithTarget);

  return (
    <main className={classes.root}>
      {Object.entries(categorizedMethodsWithTarget).map(
        ([noSteps, methodsWithTarget], index) => {
          return (
            <MethodStepAccordion
              key={noSteps}
              noSteps={Number(noSteps)}
              open={index === 0}
              methodsWithTarget={methodsWithTarget}
            />
          );
        }
      )}
    </main>
  );
};
