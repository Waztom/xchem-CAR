import React, { useState } from 'react';

import { Form } from 'react-bootstrap';
import InputGroup from 'react-bootstrap/InputGroup';
import { patchChange } from '../Utils';

const SetLayer = ({ action, updateAction }) => {
  const layer = action.layer.capitalize();
  const actiontype = action.actiontype;
  const id = action.id;

  const [Layer, setLayer] = useState(layer);

  const handleLayerChange = (e) => {
    const newlayer = e.target.value.toLowerCase();
    setLayer(e.target.value);
    patchChange(actiontype, id, 'layer', newlayer);
    updateAction(id, 'layer', newlayer);
  };

  return (
    <InputGroup size="sm" className="mb-3" key={id.toString()}>
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">Layer</InputGroup.Text>
      </InputGroup.Prepend>
      <Form.Control
        as="select"
        onChange={(event) => handleLayerChange(event)}
        size="sm"
        type="text"
        value={Layer}
      >
        <option>Organic</option>
        <option>Aqueous</option>
      </Form.Control>
    </InputGroup>
  );
};

export default SetLayer;
