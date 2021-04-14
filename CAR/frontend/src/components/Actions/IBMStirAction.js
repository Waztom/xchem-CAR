import React from "react";
import Container from "react-bootstrap/Container";
import Row from "react-bootstrap/Row";
import Col from "react-bootstrap/Col";

import SetDuration from "../SetActionInputs/SetDuration";
import SetStirring from "../SetActionInputs/SetStirring";
import SetAtmosphere from "../SetActionInputs/SetAtmosphere";
import SetTemperature from "../SetActionInputs/SetTemperature";
import IBMSetTemperatureAction from "./IBMSetTemperatureAction";

const IBMStirAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container key={actionno}>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetDuration
            action={action}
            updateAction={updateAction}
          ></SetDuration>
          <SetStirring
            action={action}
            updateAction={updateAction}
          ></SetStirring>
          <SetAtmosphere
            action={action}
            updateAction={updateAction}
          ></SetAtmosphere>
          <SetTemperature
            action={action}
            updateAction={updateAction}
          ></SetTemperature>
        </Col>
      </Row>
    </Container>
  );
};

export default IBMStirAction;
