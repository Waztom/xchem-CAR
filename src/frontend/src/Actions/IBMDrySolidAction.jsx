import React from 'react';
import Container from 'react-bootstrap/Container';
import Row from 'react-bootstrap/Row';
import Col from 'react-bootstrap/Col';

import SetTemperature from '../SetActionInputs/SetTemperature';
import SetDuration from '../SetActionInputs/SetDuration';
import SetAtmosphere from '../SetActionInputs/SetAtmosphere';

const IBMDrySolidAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetTemperature
            action={action}
            updateAction={updateAction}
          ></SetTemperature>
          <SetDuration
            action={action}
            updateAction={updateAction}
          ></SetDuration>
          <SetAtmosphere
            action={action}
            updateAction={updateAction}
          ></SetAtmosphere>
        </Col>
      </Row>
    </Container>
  );
};

export default IBMDrySolidAction;
