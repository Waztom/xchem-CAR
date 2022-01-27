import React from 'react';
import Container from 'react-bootstrap/Container';
import Row from 'react-bootstrap/Row';
import Col from 'react-bootstrap/Col';

import SetMaterial from '../SetActionInputs/SetMaterial';
import SetQuantity from '../SetActionInputs/SetQuantity.jsx';
import SetDropwise from '../SetActionInputs/SetDropwise';
import SetTemperature from '../SetActionInputs/SetTemperature';

const IBMQuenchAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container key={actionno}>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetMaterial
            action={action}
            updateAction={updateAction}
          ></SetMaterial>
          <SetQuantity
            action={action}
            updateAction={updateAction}
            name={'material'}
          ></SetQuantity>
          <SetDropwise
            action={action}
            updateAction={updateAction}
          ></SetDropwise>
          <SetTemperature
            action={action}
            updateAction={updateAction}
          ></SetTemperature>
        </Col>
      </Row>
    </Container>
  );
};

export default IBMQuenchAction;
