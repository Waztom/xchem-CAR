import React, { useState } from "react";
import axios from "axios";

import { Form } from "react-bootstrap";
import InputGroup from "react-bootstrap/InputGroup";
import FormControl from "react-bootstrap/FormControl";

const SetQuantity = ({ action, updateAction, name }) => {
  const quantname = name + "quantity";
  const unitname = name + "quantityunit";

  const quantity = action[quantname];
  const unit = action[unitname];
  const actiontype = action.actiontype;
  const id = action.id;

  const [Quantity, setQuantity] = useState({ quantity });
  const [Unit, setUnit] = useState({ unit });

  async function patchQuantity(value) {
    try {
      const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
        [quantname]: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  async function patchQuantityUnit(value) {
    try {
      const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
        [unitname]: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  const handleQuantityChange = (e) => {
    const inputQuantity = e.target.value;

    if (!isNaN(inputQuantity)) {
      setQuantity(e.target.value);
      patchQuantity(Number(e.target.value));
      updateAction(id, quantname, Number(inputQuantity));
    } else {
      alert("Please input an integer value");
    }
  };

  const handleUnitChange = (e) => {
    const inputUnit = e.target.value;

    setUnit(inputUnit);
    patchQuantityUnit(inputUnit);
    updateAction(id, unitname, inputUnit);
  };

  return (
    <InputGroup size="sm" className="mb-3" key={id.toString()}>
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">Quantity</InputGroup.Text>
      </InputGroup.Prepend>
      <FormControl
        aria-label="Small"
        aria-describedby="inputGroup-sizing-sm"
        placeholder={Quantity.quantity}
        onChange={(event) => handleQuantityChange(event)}
      />
      <Form.Control
        as="select"
        onChange={(event) => handleUnitChange(event)}
        size="sm"
        type="text"
        value={Unit.unit}
      >
        <option>moleq</option>
        <option>ml</option>
        <option>mmol</option>
      </Form.Control>
    </InputGroup>
  );
};

export default SetQuantity;
