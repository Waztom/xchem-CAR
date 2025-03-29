import { Button, Dialog, DialogActions, DialogContent, DialogTitle } from '@material-ui/core';
import React from 'react';

export const SubmitDialog = ({
  open,
  onClose,
  onCancel,
  cancelDisabled,
  CancelButtonProps = {},
  onSubmit,
  submitDisabled,
  SubmitButtonProps = {},
  id,
  title,
  content,
  ...other
}) => {
  return (
    <Dialog aria-labelledby={id} open={open} onClose={onClose} fullWidth {...other}>
      <DialogTitle id={id}>{title}</DialogTitle>
      <DialogContent dividers>{content}</DialogContent>
      <DialogActions>
        <Button
          autoFocus
          onClick={onCancel || onClose}
          color="primary"
          disabled={cancelDisabled}
          children="
          Cancel"
          {...CancelButtonProps}
        />
        <Button onClick={onSubmit} color="primary" disabled={submitDisabled} children="Submit" {...SubmitButtonProps} />
      </DialogActions>
    </Dialog>
  );
};
