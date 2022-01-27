import React from 'react';
import Overlay from 'react-bootstrap/Overlay';
import Tooltip from 'react-bootstrap/Tooltip';

function IntegerWarning({ show, target, hideTooltip, placement }) {
  return (
    <>
      <Overlay
        target={target.current}
        show={show}
        placement={placement}
        onEntered={() => hideTooltip()}
      >
        <Tooltip id="overlay-example">Enter a postive integer value</Tooltip>
      </Overlay>
    </>
  );
}

export default IntegerWarning;
