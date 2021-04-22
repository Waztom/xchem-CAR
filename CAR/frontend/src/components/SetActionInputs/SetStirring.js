import React, { useState } from "react";
import axios from "axios";

import { Form } from "react-bootstrap";
import InputGroup from "react-bootstrap/InputGroup";

const SetStirring = ({ action, updateAction }) => {
  const stirringspeed = action.stirringspeed;
  const actiontype = action.actiontype;
  const id = action.id;

  const [StirringSpeed, setStirringSpeed] = useState(stirringspeed);

  async function patchStirringSpeed(value) {
    try {
      const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
        stirringspeed: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  const handleStirringSpeedChange = (e) => {
    const newspeed = e.target.value;
    setStirringSpeed(newspeed);
    patchStirringSpeed(nwspeed);
    updateAction(id, "stirringspeed", newspeed);
  };

  return (
    <InputGroup size="sm" className="mb-3" key={id.toString()}>
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">
          Stirring speed
        </InputGroup.Text>
      </InputGroup.Prepend>
      <Form.Control
        as="select"
        onChange={(event) => handleStirringSpeedChange(event)}
        size="sm"
        type="text"
        value={StirringSpeed}
      >
        <option>gentle</option>
        <option>normal</option>
        <option>vigourous</option>
      </Form.Control>
    </InputGroup>
  );
};

export default SetStirring;
