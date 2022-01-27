import React, { useState, useRef } from 'react';

import InputGroup from 'react-bootstrap/InputGroup';
import FormControl from 'react-bootstrap/FormControl';
import IntegerWarning from '../TooltipsWarnings/IntegerWarning';
import { isInt, patchChange } from '../Utils';

const SetNumberRepetitions = ({ action, updateAction }) => {
  const numberrepetitions = action.numberofrepetitions;
  const actiontype = action.actiontype;
  const id = action.id;
  const target = useRef(null);

  const [NumberRepetitions, setNumberRepetitions] = useState(numberrepetitions);
  const [Show, setShow] = useState(false);

  function hideTooltip() {
    setTimeout(() => setShow(false), 2000);
  }

  const handleNumberRepetitionsChange = (e) => {
    const inputQuantity = Number(e.target.value);

    if (isInt(inputQuantity)) {
      setNumberRepetitions(inputQuantity);
      patchChange(actiontype, id, 'numberofrepetitions', inputQuantity);
      updateAction(id, 'numberofrepetitions', inputQuantity);
    } else {
      setNumberRepetitions(NumberRepetitions);
      setShow(true);
    }
  };

  return (
    <InputGroup size="sm" className="mb-3" key={id.toString()}>
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">
          No repetitions
        </InputGroup.Text>
      </InputGroup.Prepend>
      <FormControl
        aria-label="Small"
        aria-describedby="inputGroup-sizing-sm"
        value={NumberRepetitions}
        onChange={(event) => handleNumberRepetitionsChange(event)}
        ref={target}
      />
      <IntegerWarning
        show={Show}
        target={target}
        hideTooltip={hideTooltip}
      ></IntegerWarning>
    </InputGroup>
  );
};

export default SetNumberRepetitions;
