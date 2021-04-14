import React, { useState } from "react";
import axios from "axios";

import InputGroup from "react-bootstrap/InputGroup";
import FormControl from "react-bootstrap/FormControl";

const SetNumberRepetitions = ({ action, updateAction }) => {
  const numberrepetitions = action.numberofrepetitions;
  const actiontype = action.actiontype;
  const id = action.id;

  const [NumberRepetitions, setNumberRepetitions] = useState({
    numberrepetitions,
  });

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

    if (!isNaN(inputQuantity)) {
      setNumberRepetitions(inputQuantity);
      patchNumberRepetitions(Number(inputQuantity));
      updateAction(id, "numberofrepetitions", inputQuantity);
    } else {
      alert("Please input an integer value");
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
        placeholder={NumberRepetitions.numberrepetitions}
        onChange={(event) => handleNumberRepetitionsChange(event)}
      />
    </InputGroup>
  );
};

export default SetNumberRepetitions;
