import React, { useState, useRef } from 'react';

import InputGroup from 'react-bootstrap/InputGroup';
import FormControl from 'react-bootstrap/FormControl';
import IntegerWarning from '../TooltipsWarnings/IntegerWarning';
import { isInt, patchChange } from '../Utils';

const SetTemperature = ({ action, updateAction }) => {
  const temperature = action.temperature;
  const actiontype = action.actiontype;
  const id = action.id;
  const target = useRef(null);

  const [Temperature, setTemperature] = useState(temperature);
  const [Show, setShow] = useState(false);

  function hideTooltip() {
    setTimeout(() => setShow(false), 2000);
  }

  const handleTemperatureChange = (e) => {
    const inputTemperature = Number(e.target.value);

    if (isInt(inputTemperature)) {
      setTemperature(inputTemperature);
      patchChange(actiontype, id, 'temperature', inputTemperature);
      updateAction(id, 'temperature', inputTemperature);
    } else {
      setTemperature(Temperature);
      setShow(true);
    }
  };

  return (
    <InputGroup size="sm" className="mb-3" key={id.toString()}>
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">Temperature</InputGroup.Text>
      </InputGroup.Prepend>
      <FormControl
        aria-label="Small"
        aria-describedby="inputGroup-sizing-sm"
        value={Temperature}
        onChange={(event) => handleTemperatureChange(event)}
        ref={target}
      />
      <IntegerWarning
        show={Show}
        target={target}
        hideTooltip={hideTooltip}
        placement="top"
      ></IntegerWarning>
    </InputGroup>
  );
};

export default SetTemperature;
