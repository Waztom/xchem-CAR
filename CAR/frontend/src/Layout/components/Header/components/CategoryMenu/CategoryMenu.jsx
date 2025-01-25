import React from 'react';
import { Menu } from '@material-ui/core';
import { SuspenseWithBoundary } from '../../../../../common/components/SuspenseWithBoundary';

export const CategoryMenu = ({ id, anchorEl, onClose, children }) => {
  return (
    <Menu
      id={id}
      anchorEl={anchorEl}
      open={!!anchorEl}
      onClose={onClose}
      anchorOrigin={{ vertical: 'bottom', horizontal: 'left' }}
      getContentAnchorEl={null}
    >
      <SuspenseWithBoundary>{children}</SuspenseWithBoundary>
    </Menu>
  );
};
