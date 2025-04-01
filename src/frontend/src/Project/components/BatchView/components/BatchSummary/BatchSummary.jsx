import React from 'react';
import { Tooltip, Typography } from '@mui/material';
import { styled } from '@mui/material/styles';
import { useGetBatchSummary } from './hooks/useGetBatchSummary';
import { IconComponent } from '../../../../../common/components/IconComponent';
import { FaFlask } from 'react-icons/fa';
import { FindInPage } from '@mui/icons-material';
import { SuspenseWithBoundary } from '../../../../../common/components/SuspenseWithBoundary';

const SummaryContainer = styled('div')(({ theme }) => ({
  display: 'flex',
  gap: theme.spacing(2)
}));

const CategoryInfo = styled('div')(({ theme }) => ({
  display: 'grid',
  gridTemplateColumns: 'repeat(2, auto)',
  alignItems: 'center',
  gap: theme.spacing(1/2)
}));

const BatchSummaryContent = () => {
  const { targets, methods } = useGetBatchSummary();

  return (
    <SummaryContainer>
      <Tooltip title={`There are ${targets} targets in total`}>
        <CategoryInfo>
          <Typography>{targets}</Typography>
          <IconComponent Component={FindInPage} />
        </CategoryInfo>
      </Tooltip>

      <Tooltip title={`There are ${methods} methods in total`}>
        <CategoryInfo>
          <Typography>{methods}</Typography>
          <IconComponent Component={FaFlask} />
        </CategoryInfo>
      </Tooltip>
    </SummaryContainer>
  );
};

export const BatchSummary = () => (
  <SuspenseWithBoundary>
    <BatchSummaryContent />
  </SuspenseWithBoundary>
);

BatchSummary.displayName = 'BatchSummary';
