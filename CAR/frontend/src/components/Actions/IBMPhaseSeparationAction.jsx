import React from 'react';
import Container from 'react-bootstrap/Container';

const IBMPhaseSeparationAction = ({ action, actionno }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container>
      <h5>
        {actionno}. {actiontype}
      </h5>
    </Container>
  );
};

export default IBMPhaseSeparationAction;
