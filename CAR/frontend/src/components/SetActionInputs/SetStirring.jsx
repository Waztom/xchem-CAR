import React, { useState } from 'react';

import { Form } from 'react-bootstrap';
import InputGroup from 'react-bootstrap/InputGroup';
import { patchChange } from '../Utils';

const SetStirring = ({ action, updateAction }) => {
  const stirringspeed = action.stirringspeed;
  const actiontype = action.actiontype;
  const id = action.id;

  const [StirringSpeed, setStirringSpeed] = useState(stirringspeed);

  const handleStirringSpeedChange = (e) => {
    const newspeed = e.target.value;
    setStirringSpeed(newspeed);
    patchChange(actiontype, id, 'stirringspeed', newspeed);
    updateAction(id, 'stirringspeed', newspeed);
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
        <option>vigorous</option>
      </Form.Control>
    </InputGroup>
  );
};

export default SetStirring;
