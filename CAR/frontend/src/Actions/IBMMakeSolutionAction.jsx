import React from 'react';
import Container from 'react-bootstrap/Container';
import Row from 'react-bootstrap/Row';
import Col from 'react-bootstrap/Col';
import Image from 'react-bootstrap/Image';

import SetQuantity from '../SetActionInputs/SetQuantity.jsx';

const IBMMakeSolutionAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <Image
            width={350}
            height={350}
            src={action.soluteimage}
            alt={action.solute}
            fluid
          />
          <SetQuantity
            action={action}
            updateAction={updateAction}
            name={'solute'}
          ></SetQuantity>
          <Image
            width={350}
            height={350}
            src={action.solventimage}
            alt={action.solvent}
            fluid
          />
          <SetQuantity
            action={action}
            updateAction={updateAction}
            name={'solvent'}
          ></SetQuantity>
        </Col>
      </Row>
    </Container>
  );
};

export default IBMMakeSolutionAction;
