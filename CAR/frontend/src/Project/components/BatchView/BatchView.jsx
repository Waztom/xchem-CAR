import React from 'react';
import { useBatchContext } from '../../hooks/useBatchContext';
import { SuspenseWithBoundary } from '../../../common/components/SuspenseWithBoundary';
import { ContentBox } from '../../../common/components/ContentBox';
import { BatchSummary } from './components/BatchSummary';
import { makeStyles } from '@material-ui/core';
import { setBatchViewRef } from '../../../common/stores/batchViewsRefsStore';
import { BatchTable } from '../BatchTable';

const useStyles = makeStyles(theme => ({
  box: {
    scrollMarginTop: `${theme.spacing(2)}px`,
    // Title wrapper
    '& > :first-child': {
      position: 'sticky',
      top: 0,
      zIndex: 100
    }
  }
}));

export const BatchView = () => {
  const classes = useStyles();

  const batch = useBatchContext();

  return (
    <ContentBox
      ref={element => setBatchViewRef(batch.id, element)}
      title={batch.batchtag}
      endAdornment={
        <SuspenseWithBoundary
          SuspenseProps={{ fallback: null }}
          ErrorBoundaryProps={{ fallbackRender: () => null }}
          enableLegacySuspense
        >
          <BatchSummary />
        </SuspenseWithBoundary>
      }
      PaperProps={{ className: classes.box }}
    >
      <BatchTable />
    </ContentBox>
  );
};
