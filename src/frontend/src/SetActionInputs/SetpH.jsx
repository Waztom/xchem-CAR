import React, { useState, useRef } from 'react';

import InputGroup from 'react-bootstrap/InputGroup';
import FormControl from 'react-bootstrap/FormControl';
import IntegerWarning from '../TooltipsWarnings/IntegerWarning';
import { isFloat, isInt, patchChange } from '../Utils';

const SetpH = ({ action, updateAction }) => {
  const ph = action.pH;
  const actiontype = action.actiontype;
  const id = action.id;
  const target = useRef(null);

  const [pH, setpH] = useState(ph);
  const [Show, setShow] = useState(false);

  function hideTooltip() {
    setTimeout(() => setShow(false), 2000);
  }

  const handlepHChange = (e) => {
    const inputQuantity = e.target.value;

    if (isFloat(Number(inputQuantity)) || isInt(Number(inputQuantity))) {
      setpH(inputQuantity);
      patchChange(actiontype, id, 'pH', inputQuantity);
      updateAction(id, 'pH', inputQuantity);
    } else {
      setpH(pH);
      setShow(true);
    }
  };

  return (
    <InputGroup size="sm" className="mb-3" key={id.toString()}>
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">pH</InputGroup.Text>
      </InputGroup.Prepend>
      <FormControl
        aria-label="Small"
        aria-describedby="inputGroup-sizing-sm"
        placeholder={pH}
        onChange={(event) => handlepHChange(event)}
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

export default SetpH;
