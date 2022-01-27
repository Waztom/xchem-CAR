import React from 'react';
import Container from 'react-bootstrap/Container';
import Row from 'react-bootstrap/Row';
import Col from 'react-bootstrap/Col';

import SetSolvent from '../SetActionInputs/SetSolvent';
import SetQuantity from '../SetActionInputs/SetQuantity.jsx';
import SetNumberRepetitions from '../SetActionInputs/SetNumberRepetitions';

const IBMExtractAction = ({ action, actionno, updateAction }) => {
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
            name={''}
          ></SetSolvent>
          <SetQuantity
            action={action}
            updateAction={updateAction}
            name={'solvent'}
          ></SetQuantity>
          <SetNumberRepetitions
            action={action}
            updateAction={updateAction}
          ></SetNumberRepetitions>
        </Col>
      </Row>
    </Container>
  );
};

export default IBMExtractAction;
