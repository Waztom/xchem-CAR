import React, { useState } from "react";
import axios from "axios";

import { Form } from "react-bootstrap";
import InputGroup from "react-bootstrap/InputGroup";

const SetAtmosphere = ({ action, updateAction }) => {
  const atmosphere = action.atmosphere.capitalize();
  const actiontype = action.actiontype;
  const id = action.id;

  const [Atmosphere, setAtmosphere] = useState({ atmosphere });

  async function patchAtmosphere(value) {
    try {
      const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
        atmosphere: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  const handleAtmosphereChange = (e) => {
    const newatmosphere = e.target.value.toLowerCase();
    setAtmosphere(e.target.value);
    patchAtmosphere(newatmosphere);
    updateAction(id, "atmosphere", newatmosphere);
  };

  return (
    <InputGroup size="sm" className="mb-3" key={id.toString()}>
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">Atmosphere</InputGroup.Text>
      </InputGroup.Prepend>
      <Form.Control
        as="select"
        onChange={(event) => handleAtmosphereChange(event)}
        size="sm"
        type="text"
        value={Atmosphere.atmosphere}
      >
        <option>Nitrogen</option>
        <option>Air</option>
      </Form.Control>
    </InputGroup>
  );
};

export default SetAtmosphere;
