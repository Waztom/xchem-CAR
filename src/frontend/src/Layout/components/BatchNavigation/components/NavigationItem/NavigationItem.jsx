import React, { useCallback } from 'react';
import { TreeItem } from '@material-ui/lab';
import { Checkbox, CircularProgress, Fab, makeStyles, Tooltip, Typography } from '@material-ui/core';
import { setBatchSelected, useBatchNavigationStore } from '../../../../../common/stores/batchNavigationStore';
import { useBatchViewsRefs } from '../../../../../common/stores/batchViewsRefsStore';
import { CgArrowsScrollV } from 'react-icons/cg';
import { IconComponent } from '../../../../../common/components/IconComponent';
import classNames from 'classnames';
import { useTemporaryId } from '../../../../../common/hooks/useTemporaryId';
import { DeleteForever } from '@material-ui/icons';
import { requestDeleteSubBatch } from '../../../../stores/deleteSubBatchDialogStore';

const useStyles = makeStyles(theme => ({
  label: {
    display: 'flex',
    minWidth: 0,
    alignItems: 'center'
  },
  action: {
    padding: 0
  },
  leaf: {
    cursor: 'default'
  },
  name: {
    flexGrow: 1
  },
  button: {
    minHeight: 'unset',
    width: theme.spacing(2),
    height: theme.spacing(2),
    boxShadow: 'none !important'
  },
  icon: {
    width: '1.2em !important',
    height: '1.2em !important'
  },
  actions: {
    display: 'flex',
    gap: theme.spacing(1 / 4),
    alignItems: 'center'
  },
  deleteButton: {
    color: theme.palette.error.main
  }
}));

export const NavigationItem = ({ node, children }) => {
  const classes = useStyles();

  const { batch, children: subBatchNodes } = node;

  const displayed = useBatchNavigationStore(useCallback(state => state.selected[batch.id] || false, [batch.id]));
  const elementRef = useBatchViewsRefs(useCallback(state => state.refs[batch.id], [batch.id]));

  // Used when newly created batch is loading
  const isTemporaryBatch = useTemporaryId().isTemporaryId(batch.id);

  // Only subBatches without children subBatches can be deleted
  const deleteEnabled = !!batch.batch_id && !subBatchNodes.length;

  return (
    <TreeItem
      classes={{ label: classes.label, content: !children.length && classes.leaf }}
      nodeId={String(batch.id)}
      label={
        <>
          <Typography className={classes.name} noWrap>
            {batch.batchtag}
          </Typography>
          <div className={classes.actions}>
            {isTemporaryBatch ? (
              <CircularProgress className={classes.icon} />
            ) : (
              <>
                {deleteEnabled && (
                  <Tooltip title="Delete batch">
                    <Fab
                      className={classNames(classes.action, classes.button, classes.deleteButton)}
                      size="small"
                      onClick={() => requestDeleteSubBatch(batch)}
                    >
                      <DeleteForever className={classes.icon} fontSize="inherit" />
                    </Fab>
                  </Tooltip>
                )}
                {!!elementRef && (
                  <Tooltip title="Scroll to batch">
                    <Fab
                      className={classNames(classes.action, classes.button)}
                      size="small"
                      onClick={event => {
                        event.stopPropagation();
                        elementRef.scrollIntoView({ behavior: 'smooth' });
                      }}
                      color="secondary"
                    >
                      <IconComponent className={classes.icon} Component={CgArrowsScrollV} />
                    </Fab>
                  </Tooltip>
                )}
                <Tooltip title={displayed ? 'Hide batch' : 'Display batch'}>
                  <Checkbox
                    checked={displayed}
                    className={classes.action}
                    onClick={e => e.stopPropagation()}
                    onChange={(_, checked) => setBatchSelected(batch.id, checked)}
                    inputProps={{ 'aria-label': batch.batchtag }}
                  />
                </Tooltip>
              </>
            )}
          </div>
        </>
      }
    >
      {children}
    </TreeItem>
  );
};
