import React from 'react';
import Container from 'react-bootstrap/Container';
import Row from 'react-bootstrap/Row';
import Col from 'react-bootstrap/Col';

import SetMaterial from '../SetActionInputs/SetMaterial';
import SetQuantity from '../SetActionInputs/SetQuantity.jsx';
import SetNumberRepetitions from '../SetActionInputs/SetNumberRepetitions';

const IBMWashAction = ({ action, actionno, updateAction }) => {
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
          <SetNumberRepetitions
            action={action}
            updateAction={updateAction}
          ></SetNumberRepetitions>
        </Col>
      </Row>
    </Container>
  );
};

export default IBMWashAction;
