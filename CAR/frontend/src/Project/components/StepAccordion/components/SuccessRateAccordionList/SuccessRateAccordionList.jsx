import { Fragment } from 'react';
import { Divider, List, ListItem, makeStyles } from '@material-ui/core';
import { useGetMethodsReactions } from './hooks/useGetMethodsReactions';
import { useCategorizeMethodsDataBySuccessRate } from './hooks/useCategorizeMethodsDataBySuccessRate';
import { LoadingSpinner } from '../../../../../common/components/LoadingSpinner';
import { SuccessRateAccordion } from '../../../SuccessRateAccordion';

const useStyles = makeStyles((theme) => ({
  list: {
    width: '100%',
    '& > li': {
      padding: 0,
    },
  },
}));

// Separated from StepAccordion to enable loading reactions only when StepAccordion's details are open
export const SuccessRateAccordionList = ({ noSteps, methodsWithTarget }) => {
  const classes = useStyles();

  const { methodsData, isLoading } = useGetMethodsReactions(methodsWithTarget);
  const categorizedMethodsData =
    useCategorizeMethodsDataBySuccessRate(methodsData);

  if (isLoading) {
    return <LoadingSpinner />;
  }

  return (
    <List className={classes.list} disablePadding>
      {Object.entries(categorizedMethodsData)
        .sort(([keyA], [keyB]) => keyB.localeCompare(keyA))
        .map(([successString, methodData], index) => {
          return (
            <Fragment key={successString}>
              <ListItem disableGutters>
                <SuccessRateAccordion
                  noSteps={noSteps}
                  successString={successString}
                  methodData={methodData}
                />
              </ListItem>
              {!!(index < methodData.length - 1) && <Divider />}
            </Fragment>
          );
        })}
    </List>
  );
};
