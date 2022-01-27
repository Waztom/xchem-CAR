import React from 'react';
import Container from 'react-bootstrap/Container';
import Row from 'react-bootstrap/Row';
import Col from 'react-bootstrap/Col';

import SetDryingAgent from '../SetActionInputs/SetDryingAgent';

const IBMDrySolutionAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetDryingAgent
            action={action}
            updateAction={updateAction}
          ></SetDryingAgent>
        </Col>
      </Row>
    </Container>
  );
};

export default IBMDrySolutionAction;
