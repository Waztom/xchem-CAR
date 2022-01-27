import React from 'react';
import Col from 'react-bootstrap/Col';

import SetGas from '../SetActionInputs/SetGas';
import SetDuration from '../SetActionInputs/SetDuration';

const IBMDegasAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetGas action={action} updateAction={updateAction}></SetGas>
          <SetDuration
            action={action}
            updateAction={updateAction}
          ></SetDuration>
        </Col>
      </Row>
    </Container>
  );
};

export default IBMDegasAction;
