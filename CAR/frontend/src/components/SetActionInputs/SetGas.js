import React, { useState } from "react";
import axios from "axios";

import InputGroup from "react-bootstrap/InputGroup";
import FormControl from "react-bootstrap/FormControl";

const SetGas = ({ action, updateAction }) => {
  const gas = action.gas;
  const actiontype = action.actiontype;
  const id = action.id;

  const [Gas, setGas] = useState(gas);

  async function patchGas(value) {
    try {
      const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
        gas: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  const handleGasChange = (e) => {
    const newgas = e.target.value;
    setGas(newgas);
    patchGas(newgas);
    updateAction(id, "gas", newgas);
  };

  return (
    <InputGroup size="sm" className="mb-3" key={id.toString()}>
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">Quantity</InputGroup.Text>
      </InputGroup.Prepend>
      <FormControl
        aria-label="Small"
        aria-describedby="inputGroup-sizing-sm"
        placeholder={Gas}
        onChange={(event) => handleGasChange(event)}
      />
    </InputGroup>
  );
};

export default SetGas;
