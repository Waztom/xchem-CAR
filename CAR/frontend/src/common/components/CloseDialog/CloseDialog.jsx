import React from 'react';
import { Button, Dialog, DialogActions, DialogContent, DialogTitle } from '@material-ui/core';

export const CloseDialog = ({ open, onClose, closeDisabled, CloseButtonProps = {}, id, title, content, ...other }) => {
  return (
    <Dialog aria-labelledby={id} open={open} onClose={onClose} {...other}>
      <DialogTitle id={id}>{title}</DialogTitle>
      <DialogContent dividers>{content}</DialogContent>
      <DialogActions>
        <Button
          autoFocus
          onClick={onClose}
          color="primary"
          disabled={closeDisabled}
          children="Close"
          {...CloseButtonProps}
        />
      </DialogActions>
    </Dialog>
  );
};
