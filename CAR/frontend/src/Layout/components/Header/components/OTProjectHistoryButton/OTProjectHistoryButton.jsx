import React, { useState } from 'react';
import { Button } from '@material-ui/core';
import { CategoryMenu } from '../CategoryMenu';
import { OTProjectMenuContents } from '../OTProjectMenuContents';
import { requestOtProjectSummary } from '../../../../stores/otProjectSummaryDialogStore';
import { OTProjectSummaryDialog } from '../../../OTProjectSummaryDialog';

export const OTProjectHistoryButton = () => {
  const [menuAnchorEl, setMenuAnchorEl] = useState(null);

  return (
    <>
      <Button onClick={event => setMenuAnchorEl(event.currentTarget)}>OT protocol history</Button>
      <CategoryMenu id="ot-protocol-history-menu" anchorEl={menuAnchorEl} onClose={() => setMenuAnchorEl(null)}>
        <OTProjectMenuContents
          onSelected={otProject => {
            setMenuAnchorEl(null);
            requestOtProjectSummary(otProject.id);
          }}
        />
      </CategoryMenu>

      <OTProjectSummaryDialog />
    </>
  );
};
