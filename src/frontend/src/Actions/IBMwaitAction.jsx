import React from 'react';
import Container from 'react-bootstrap/Container';
import Row from 'react-bootstrap/Row';
import Col from 'react-bootstrap/Col';

import SetDuration from '../SetActionInputs/SetDuration';
import SetTemperature from '../SetActionInputs/SetTemperature';

const IBMWaitAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container key={actionno}>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetDuration
            action={action}
            updateAction={updateAction}
          ></SetDuration>
          <SetTemperature
            action={action}
            updateAction={updateAction}
          ></SetTemperature>
        </Col>
      </Row>
    </Container>
  );
};

export default IBMWaitAction;
