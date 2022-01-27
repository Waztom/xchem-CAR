import React from 'react';
import Container from 'react-bootstrap/Container';
import Row from 'react-bootstrap/Row';
import Col from 'react-bootstrap/Col';
import Image from 'react-bootstrap/Image';

const IBMStoreAction = ({ action, actionno }) => {
  const actiontype = action.actiontype.capitalize();
  const material = action.material.capitalize();

  return (
    <Container key={actionno}>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <Image src={action.materialimage} alt={action.material} fluid />
        </Col>
      </Row>
    </Container>
  );
};

export default IBMStoreAction;
