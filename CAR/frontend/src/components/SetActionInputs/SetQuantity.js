import React, { useState, useRef } from "react";
import axios from "axios";

import { Form } from "react-bootstrap";
import InputGroup from "react-bootstrap/InputGroup";
import IntegerWarning from "../Tooltips/IntegerWarning";

const SetQuantity = ({ action, updateAction, name }) => {
  const quantname = name + "quantity";
  const unitname = name + "quantityunit";
  const target = useRef(null);

  const quantity = action[quantname];
  const unit = action[unitname];
  const actiontype = action.actiontype;
  const id = action.id;

  const [Quantity, setQuantity] = useState(quantity);
  const [Unit, setUnit] = useState(unit);
  const [Show, setShow] = useState(false);

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
      setQuantity(inputQuantity);
      patchQuantity(Number(inputQuantity));
      updateAction(id, quantname, Number(inputQuantity));
    } else {
      setQuantity(Quantity);
      setShow(true);
    }
  };

  // Could not find a way to handle timers better - without
  // passing this function to child component...
  function hideTooltip() {
    setTimeout(() => setShow(false), 2000);
  }

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
      <Form.Control
        aria-label="Small"
        aria-describedby="inputGroup-sizing-sm"
        value={Quantity}
        onChange={(event) => handleQuantityChange(event)}
        ref={target}
      />
      <IntegerWarning
        show={Show}
        target={target}
        hideTooltip={hideTooltip}
        placement="top"
      ></IntegerWarning>
      <Form.Control
        as="select"
        onChange={(event) => handleUnitChange(event)}
        size="sm"
        type="text"
        value={Unit}
      >
        <option>moleq</option>
        <option>ml</option>
        <option>mmol</option>
      </Form.Control>
    </InputGroup>
  );
};

export default SetQuantity;
