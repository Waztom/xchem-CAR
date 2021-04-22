import React, { useState, useRef } from "react";
import axios from "axios";

import InputGroup from "react-bootstrap/InputGroup";
import FormControl from "react-bootstrap/FormControl";
import IntegerWarning from "../Tooltips/IntegerWarning";

const SetNumberRepetitions = ({ action, updateAction }) => {
  const numberrepetitions = action.numberofrepetitions;
  const actiontype = action.actiontype;
  const id = action.id;
  const target = useRef(null);

  const [NumberRepetitions, setNumberRepetitions] = useState(numberrepetitions);
  const [Show, setShow] = useState(false);

  async function patchNumberRepetitions(value) {
    try {
      const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
        numberofrepetitions: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  const handleNumberRepetitionsChange = (e) => {
    const inputQuantity = e.target.value;

    if (!isNaN(inputQuantity) && Number.isInteger(inputQuantity)) {
      setNumberRepetitions(inputQuantity);
      patchNumberRepetitions(Number(inputQuantity));
      updateAction(id, "numberofrepetitions", Number(inputQuantity));
    } else {
      setShow(true);
    }
  };

  // Could not find a way to handle timers better - without
  // passing this function to child component...
  function hideTooltip() {
    setTimeout(() => setShow(false), 2000);
  }

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
        placeholder={NumberRepetitions}
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
