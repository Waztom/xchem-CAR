import React from 'react';
import { styled } from '@mui/material/styles';
import { useBatchContext } from '../../hooks/useBatchContext';
import { SuspenseWithBoundary } from '../../../common/components/SuspenseWithBoundary';
import { ContentBox } from '../../../common/components/ContentBox';
import { BatchSummary } from './components/BatchSummary';
import { setBatchViewRef } from '../../../common/stores/batchViewsRefsStore';
import { BatchTable } from '../BatchTable';

const StyledContentBox = styled(ContentBox)(({ theme }) => ({
  scrollMarginTop: theme.spacing(2),
  '& > :first-of-type': {
    position: 'sticky',
    top: 0,
    zIndex: theme.zIndex.appBar
  }
}));

const BatchViewContent = () => {
  const batch = useBatchContext();
  
  return (
    <StyledContentBox
      ref={element => setBatchViewRef(batch.id, element)}
      title={batch.batchtag}
      endAdornment={
        <SuspenseWithBoundary
          SuspenseProps={{ fallback: null }}
          ErrorBoundaryProps={{ fallbackRender: () => null }}
        >
          <BatchSummary />
        </SuspenseWithBoundary>
      }
    >
      <SuspenseWithBoundary>
        <BatchTable />
      </SuspenseWithBoundary>
    </StyledContentBox>
  );
};

export const BatchView = () => {
  return (
    <SuspenseWithBoundary>
      <BatchViewContent />
    </SuspenseWithBoundary>
  );
};

BatchView.displayName = 'BatchView';
