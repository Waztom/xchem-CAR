import React from 'react';
import Container from 'react-bootstrap/Container';
import Row from 'react-bootstrap/Row';
import Col from 'react-bootstrap/Col';

import SetPhaseTokeep from '../SetActionInputs/SetPhaseTokeep';
import SetSolvent from '../SetActionInputs/SetSolvent';
import SetQuantity from '../SetActionInputs/SetQuantity.jsx';

const IBMFilterAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetPhaseTokeep
            action={action}
            updateAction={updateAction}
          ></SetPhaseTokeep>
          <SetSolvent
            action={action}
            updateAction={updateAction}
            name={'rinsing'}
          ></SetSolvent>
          <SetQuantity
            action={action}
            updateAction={updateAction}
            name={'rinsingsolvent'}
          ></SetQuantity>
          <SetSolvent
            action={action}
            updateAction={updateAction}
            name={'extractionforprecipitate'}
          ></SetSolvent>
          <SetQuantity
            action={action}
            updateAction={updateAction}
            name={'extractionforprecipitatesolvent'}
          ></SetQuantity>
        </Col>
      </Row>
    </Container>
  );
};

export default IBMFilterAction;
