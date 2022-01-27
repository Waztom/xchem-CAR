import React, { useState, useRef } from 'react';

import { Form } from 'react-bootstrap';
import InputGroup from 'react-bootstrap/InputGroup';
import IntegerWarning from '../TooltipsWarnings/IntegerWarning';
import { isFloat, isInt, patchChange } from '../Utils';

const SetDuration = ({ action, updateAction }) => {
  const duration = action.duration;
  const unit = action.durationunit;
  const actiontype = action.actiontype;
  const id = action.id;
  const target = useRef(null);

  const [Duration, setDuration] = useState(duration);
  const [Unit, setUnit] = useState(unit);
  const [Show, setShow] = useState(false);

  function hideTooltip() {
    setTimeout(() => setShow(false), 2000);
  }

  const handleDurationChange = (e) => {
    const inputDuration = e.target.value;

    if (isFloat(Number(inputDuration)) || isInt(Number(inputDuration))) {
      setDuration(inputDuration);
      patchChange(actiontype, id, 'duration', inputDuration);
      updateAction(id, 'duration', Number(inputDuration));
    } else {
      setDuration(Duration);
      setShow(true);
    }
  };

  const handleUnitChange = (e) => {
    const newunit = e.target.value;
    setUnit(newunit);
    patchChange(actiontype, id, 'durationunit', newunit);
    updateAction(id, 'durationunit', newunit);
  };

  return (
    <InputGroup size="sm" className="mb-3" key={id.toString()}>
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">Duration</InputGroup.Text>
      </InputGroup.Prepend>
      <Form.Control
        aria-label="Small"
        aria-describedby="inputGroup-sizing-sm"
        onChange={(event) => handleDurationChange(event)}
        value={Duration}
        ref={target}
      />
      <IntegerWarning
        show={Show}
        target={target}
        hideTooltip={hideTooltip}
        postion="right"
      ></IntegerWarning>
      <Form.Control
        as="select"
        onChange={(event) => handleUnitChange(event)}
        size="sm"
        type="text"
        value={Unit}
      >
        <option>seconds</option>
        <option>minutes</option>
        <option>hours</option>
      </Form.Control>
    </InputGroup>
  );
};

export default SetDuration;
