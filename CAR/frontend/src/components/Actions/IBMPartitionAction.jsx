import React from 'react';
import Container from 'react-bootstrap/Container';
import Row from 'react-bootstrap/Row';
import Col from 'react-bootstrap/Col';

import SetSolvent from '../SetActionInputs/SetSolvent';
import SetQuantity from '../SetActionInputs/SetQuantity.jsx';

const IBMPartitionAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetSolvent
            action={action}
            updateAction={updateAction}
            name={'firstparition'}
          ></SetSolvent>
          <SetQuantity
            action={action}
            updateAction={updateAction}
            name={'firstpartitionsolvent'}
          ></SetQuantity>
          <SetSolvent
            action={action}
            updateAction={updateAction}
            name={'secondpartition'}
          ></SetSolvent>
          <SetQuantity
            action={action}
            updateAction={updateAction}
            name={'secondpartitionsolvent'}
          ></SetQuantity>
        </Col>
      </Row>
    </Container>
  );
};

export default IBMPartitionAction;
