import React from 'react';
import Container from 'react-bootstrap/Container';
import Row from 'react-bootstrap/Row';
import Col from 'react-bootstrap/Col';

import SetTemperature from '../SetActionInputs/SetTemperature';

const IBMSetTemperatureAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container key={actionno}>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetTemperature
            action={action}
            updateAction={updateAction}
          ></SetTemperature>
        </Col>
      </Row>
    </Container>
  );
};

export default IBMSetTemperatureAction;
