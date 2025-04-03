import React, { useState } from 'react';

import { Form } from 'react-bootstrap';
import InputGroup from 'react-bootstrap/InputGroup';
import { patchChange } from '../Utils';

const SetDeanStark = ({ action, updateAction }) => {
  const deanstark = action.deanstarkapparatus;
  const actiontype = action.actiontype;
  const id = action.id;

  const [DeanStark, setDeanStark] = useState(deanstark);

  const handleDeanStarkChange = (e) => {
    const newdeanstark = e.target.value;
    setDeanStark(newdeanstark);
    patchChange(actiontype, id, 'deanstarkapparatus', newdeanstark);
    updateAction(id, 'deanstarkapparatus', newdeanstark);
  };

  return (
    <InputGroup size="sm" className="mb-3" key={id.toString()}>
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">Dean Stark</InputGroup.Text>
      </InputGroup.Prepend>
      <Form.Control
        as="select"
        onChange={(event) => handleDeanStarkChange(event)}
        size="sm"
        type="text"
        value={DeanStark}
      >
        <option>True</option>
        <option>False</option>
      </Form.Control>
    </InputGroup>
  );
};

export default SetDeanStark;
