import React from 'react';
import Container from 'react-bootstrap/Container';
import Row from 'react-bootstrap/Row';
import Col from 'react-bootstrap/Col';

import SetLayer from '../SetActionInputs/SetLayer';

const IBMCollectLayerAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetLayer action={action} updateAction={updateAction}></SetLayer>
        </Col>
      </Row>
    </Container>
  );
};

export default IBMCollectLayerAction;
