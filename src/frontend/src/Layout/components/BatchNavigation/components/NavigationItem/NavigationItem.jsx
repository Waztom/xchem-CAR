import React, { useCallback } from 'react';
import { TreeItem } from '@mui/x-tree-view/TreeItem';
import { Checkbox, CircularProgress, Fab, Tooltip, Typography } from '@mui/material';
import { styled } from '@mui/material/styles';
import { DeleteForever } from '@mui/icons-material';
import { CgArrowsScrollV } from 'react-icons/cg';
import { SuspenseWithBoundary } from '../../../../../common/components/SuspenseWithBoundary';
import { IconComponent } from '../../../../../common/components/IconComponent';
import { setBatchSelected, useBatchNavigationStore } from '../../../../../common/stores/batchNavigationStore';
import { useBatchViewsRefs } from '../../../../../common/stores/batchViewsRefsStore';
import { useTemporaryId } from '../../../../../common/hooks/useTemporaryId';
import { requestDeleteSubBatch } from '../../../../stores/deleteSubBatchDialogStore';

const StyledTreeItem = styled(TreeItem)(({ theme }) => ({
  '& .MuiTreeItem-label': {
    display: 'flex',
    minWidth: 0,
    alignItems: 'center'
  },
  '& .MuiTreeItem-content.Mui-leaf': {
    cursor: 'default'
  }
}));

const BatchName = styled(Typography)(() => ({
  flexGrow: 1
}));

const ActionsWrapper = styled('div')(({ theme }) => ({
  display: 'flex',
  gap: theme.spacing(1/4),
  alignItems: 'center'
}));

const StyledFab = styled(Fab)(({ theme }) => ({
  minHeight: 'unset',
  width: theme.spacing(3),  // Increased size slightly
  height: theme.spacing(3),
  boxShadow: 'none !important',
  padding: 0,
  '&:hover': {
    backgroundColor: theme.palette.error.light,
    '& .MuiSvgIcon-root': {
      color: theme.palette.error.contrastText
    }
  }
}));

const StyledIcon = styled('span')(() => ({
  width: '1.2em !important',
  height: '1.2em !important'
}));

const StyledDeleteIcon = styled(DeleteForever)(({ theme }) => ({
  fontSize: '1.2rem',
  color: theme.palette.error.main,
}));

const NavigationItemContent = ({ batch, subBatchNodes, elementRef }) => {
  const displayed = useBatchNavigationStore(useCallback(state => state.selected[batch.id] || false, [batch.id]));
  const { isTemporaryId } = useTemporaryId();
  const isTemporaryBatch = isTemporaryId(batch.id);
  const deleteEnabled = !!batch.batch_id && !subBatchNodes.length;

  return (
    <>
      <BatchName noWrap>{batch.batchtag}</BatchName>
      <ActionsWrapper>
        {isTemporaryBatch ? (
          <CircularProgress sx={{ width: '1.2em', height: '1.2em' }} />
        ) : (
          <>
            {deleteEnabled && (
              <Tooltip title="Delete batch">
                <StyledFab
                  size="small"
                  onClick={(e) => {
                    e.stopPropagation();
                    requestDeleteSubBatch(batch);
                  }}
                >
                  <StyledDeleteIcon />
                </StyledFab>
              </Tooltip>
            )}
            {!!elementRef && (
              <Tooltip title="Scroll to batch">
                <StyledFab
                  size="small"
                  onClick={event => {
                    event.stopPropagation();
                    elementRef.scrollIntoView({ behavior: 'smooth' });
                  }}
                  color="secondary"
                >
                  <IconComponent component={StyledIcon} Component={CgArrowsScrollV} />
                </StyledFab>
              </Tooltip>
            )}
            <Tooltip title={displayed ? 'Hide batch' : 'Display batch'}>
              <Checkbox
                checked={displayed}
                onClick={e => e.stopPropagation()}
                onChange={(_, checked) => setBatchSelected(batch.id, checked)}
                inputProps={{ 'aria-label': batch.batchtag }}
              />
            </Tooltip>
          </>
        )}
      </ActionsWrapper>
    </>
  );
};

export const NavigationItem = ({ node, children }) => {
  const { batch, children: subBatchNodes } = node;
  const elementRef = useBatchViewsRefs(useCallback(state => state.refs[batch.id], [batch.id]));

  return (
    <StyledTreeItem
      nodeId={String(batch.id)}
      label={
        <SuspenseWithBoundary>
          <NavigationItemContent 
            batch={batch}
            subBatchNodes={subBatchNodes}
            elementRef={elementRef}
          />
        </SuspenseWithBoundary>
      }
    >
      {children}
    </StyledTreeItem>
  );
};

NavigationItem.displayName = 'NavigationItem';
