import React, { useState } from "react";
import axios from "axios";

import { Form } from "react-bootstrap";
import InputGroup from "react-bootstrap/InputGroup";
import FormControl from "react-bootstrap/FormControl";

const SetDuration = ({ action, updateAction }) => {
  const duration = action.duration;
  const unit = action.durationunit;
  const actiontype = action.actiontype;
  const id = action.id;

  const [Duration, setDuration] = useState({ duration });
  const [Unit, setUnit] = useState({ unit });

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
    const inputQuantity = e.target.value;

    if (!isNaN(inputQuantity)) {
      setDuration(inputQuantity);
      patchDuration(Number(inputQuantity));
      updateAction(id, "duration", inputQuantity);
    } else {
      alert("Please input an integer value");
    }
  };

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
      <FormControl
        aria-label="Small"
        aria-describedby="inputGroup-sizing-sm"
        placeholder={Duration.duration}
        onChange={(event) => handleDurationChange(event)}
      />
      <Form.Control
        as="select"
        onChange={(event) => handleUnitChange(event)}
        size="sm"
        type="text"
        value={Unit.unit}
      >
        <option>seconds</option>
        <option>minutes</option>
        <option>hours</option>
      </Form.Control>
    </InputGroup>
  );
};

export default SetDuration;
