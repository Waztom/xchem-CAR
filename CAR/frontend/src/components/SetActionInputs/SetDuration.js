import React, { useState, useRef } from "react";
import axios from "axios";

import { Form } from "react-bootstrap";
import InputGroup from "react-bootstrap/InputGroup";
import IntegerWarning from "../Tooltips/IntegerWarning";

const SetDuration = ({ action, updateAction }) => {
  const duration = action.duration;
  const unit = action.durationunit;
  const actiontype = action.actiontype;
  const id = action.id;
  const target = useRef(null);

  const [Duration, setDuration] = useState(duration);
  const [Unit, setUnit] = useState(unit);
  const [Show, setShow] = useState(false);

  async function patchDuration(value) {
    try {
      const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
        duration: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  async function patchDurationUnit(value) {
    try {
      const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
        durationunit: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  const handleDurationChange = (e) => {
    const inputDuration = e.target.value;

    if (!isNaN(inputDuration)) {
      setDuration(inputDuration);
      patchDuration(Number(inputDuration));
      updateAction(id, "duration", Number(inputDuration));
    } else {
      setDuration(Duration);
      setShow(true);
    }
  };

  // Could not find a way to handle timers better - without
  // passing this function to child component...
  function hideTooltip() {
    setTimeout(() => setShow(false), 2000);
  }

  const handleUnitChange = (e) => {
    const newunit = e.target.value;
    setUnit(newunit);
    patchDurationUnit(newunit);
    updateAction(id, "durationunit", newunit);
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
