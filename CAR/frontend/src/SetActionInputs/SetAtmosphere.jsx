import React, { useState } from 'react';

import { Form } from 'react-bootstrap';
import InputGroup from 'react-bootstrap/InputGroup';
import { patchChange } from '../Utils';

const SetAtmosphere = ({ action, updateAction }) => {
  const atmosphere = action.atmosphere.capitalize();
  const actiontype = action.actiontype;
  const id = action.id;

  const [Atmosphere, setAtmosphere] = useState(atmosphere);

  const handleAtmosphereChange = (e) => {
    const newatmosphere = e.target.value.toLowerCase();
    setAtmosphere(e.target.value);
    patchChange(actiontype, id, 'atmosphere', newatmosphere);
    updateAction(id, 'atmosphere', newatmosphere);
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
        value={Atmosphere}
      >
        <option>Nitrogen</option>
        <option>Air</option>
      </Form.Control>
    </InputGroup>
  );
};

export default SetAtmosphere;
