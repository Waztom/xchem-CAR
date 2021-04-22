import React, { useState } from "react";
import axios from "axios";

import { Form } from "react-bootstrap";
import InputGroup from "react-bootstrap/InputGroup";

const SetPhaseTokeep = ({ action, updateAction }) => {
  const phasetokeep = action.phasetokeep.capitalize();
  const actiontype = action.actiontype;
  const id = action.id;

  const [PhaseToKeep, setPhaseToKeep] = useState(phasetokeep);

  async function patchPhaseToKeep(value) {
    try {
      const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
        phasetokeep: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  const handlePhaseToKeepChange = (e) => {
    const newphasetokeep = e.target.value.toLowerCase();
    setPhaseToKeep(newphasetokeep);
    patchPhaseToKeep(newphasetokeep);
    updateAction(id, "phasetokeep", newphasetokeep);
  };

  return (
    <InputGroup size="sm" className="mb-3" key={id.toString()}>
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">
          Phase to keep
        </InputGroup.Text>
      </InputGroup.Prepend>
      <Form.Control
        as="select"
        onChange={(event) => handlePhaseToKeepChange(event)}
        size="sm"
        type="text"
        value={PhaseToKeep}
      >
        <option>Filtrate</option>
        <option>Precipitate</option>
      </Form.Control>
    </InputGroup>
  );
};

export default SetPhaseTokeep;
