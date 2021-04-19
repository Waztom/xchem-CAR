import React from "react";
import Image from "react-bootstrap/Image";
import Container from "react-bootstrap/Container";
import Row from "react-bootstrap/Row";
import Col from "react-bootstrap/Col";
import Button from "react-bootstrap/Button";
import { Trash } from "react-bootstrap-icons";

import SetQuantity from "../SetActionInputs/SetQuantity";
import SetDropwise from "../SetActionInputs/SetDropwise";
import SetAtmosphere from "../SetActionInputs/SetAtmosphere";

import JSME from "../MolDrawer/JSME";

const IBMAddAction = ({ action, actionno, updateAction, handleDelete }) => {
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
      <Button
        key={action.id}
        onClick={() => handleDelete(action.actiontype, action.id)}
      >
        <Trash></Trash>
      </Button>
      <JSME></JSME>
    </Container>
  );
};

export default IBMAddAction;
