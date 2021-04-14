import React from "react";
import Image from "react-bootstrap/Image";
import Container from "react-bootstrap/Container";
import Row from "react-bootstrap/Row";
import Col from "react-bootstrap/Col";

import SetQuantity from "../SetActionInputs/SetQuantity";
import SetDropwise from "../SetActionInputs/SetDropwise";
import SetAtmosphere from "../SetActionInputs/SetAtmosphere";

const IBMAddAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container key={actionno}>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <Image
            width={350}
            height={350}
            src={action.materialimage}
            alt={action.material}
            fluid
          />
          <SetQuantity
            action={action}
            updateAction={updateAction}
            name={"material"}
          ></SetQuantity>
          <SetDropwise
            action={action}
            updateAction={updateAction}
          ></SetDropwise>
          <SetAtmosphere
            action={action}
            updateAction={updateAction}
          ></SetAtmosphere>
        </Col>
      </Row>
    </Container>
  );
};

export default IBMAddAction;
