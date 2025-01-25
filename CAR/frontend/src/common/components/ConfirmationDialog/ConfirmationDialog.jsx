import React from 'react';
import { Button, Dialog, DialogActions, DialogContent, DialogTitle } from '@material-ui/core';

export const ConfirmationDialog = ({
  open,
  onClose,
  onCancel,
  cancelDisabled,
  CancelButtonProps = {},
  onOk,
  okDisabled,
  OkButtonProps = {},
  id,
  title,
  content,
  ...other
}) => {
  return (
    <Dialog aria-labelledby={id} open={open} onClose={onClose} {...other}>
      <DialogTitle id={id}>{title}</DialogTitle>
      <DialogContent dividers>{content}</DialogContent>
      <DialogActions>
        <Button
          autoFocus
          onClick={onCancel || onClose}
          disabled={cancelDisabled}
          color="primary"
          children="Cancel"
          {...CancelButtonProps}
        />
        <Button onClick={onOk} disabled={okDisabled} color="primary" children="Ok" {...OkButtonProps} />
      </DialogActions>
    </Dialog>
  );
};
