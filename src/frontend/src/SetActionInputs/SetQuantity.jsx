import React, { useState, useRef } from 'react';

import { Form } from 'react-bootstrap';
import InputGroup from 'react-bootstrap/InputGroup';
import IntegerWarning from '../TooltipsWarnings/IntegerWarning';
import { isFloat, isInt, patchChange } from '../Utils';

const SetQuantity = ({ action, updateAction, name }) => {
  const quantname = name + 'quantity';
  const unitname = name + 'quantityunit';
  const target = useRef(null);

  const quantity = action[quantname];
  const unit = action[unitname];
  const actiontype = action.actiontype;
  const id = action.id;

  const [Quantity, setQuantity] = useState(quantity);
  const [Unit, setUnit] = useState(unit);
  const [Show, setShow] = useState(false);

  console.log(Unit);

  function hideTooltip() {
    setTimeout(() => setShow(false), 2000);
  }

  const handleQuantityChange = (e) => {
    const inputQuantity = e.target.value;

    if (isFloat(Number(inputQuantity)) || isInt(Number(inputQuantity))) {
      setQuantity(inputQuantity);
      patchChange(actiontype, id, quantname, inputQuantity);
      updateAction(id, quantname, inputQuantity);
    } else {
      setQuantity(Quantity);
      setShow(true);
    }
  };

  const handleUnitChange = (e) => {
    const inputUnit = e.target.value;
    setUnit(inputUnit);
    patchChange(actiontype, id, unitname, inputUnit);
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
        <option>ul</option>
        <option>mmol</option>
      </Form.Control>
    </InputGroup>
  );
};

export default SetQuantity;
