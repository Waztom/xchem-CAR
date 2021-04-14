import React, { useState } from "react";
import axios from "axios";

import InputGroup from "react-bootstrap/InputGroup";
import FormControl from "react-bootstrap/FormControl";

const SetTemperature = ({ action, updateAction }) => {
  const temperature = action.temperature;
  const actiontype = action.actiontype;
  const id = action.id;

  const [Temperature, setTemperature] = useState({ temperature });

  async function patchTemperature(value) {
    try {
      const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
        temperature: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  const handleTemperatureChange = (e) => {
    const inputQuantity = e.target.value;

    if (!isNaN(inputQuantity)) {
      setTemperature(inputQuantity);
      patchTemperature(Number(inputQuantity));
      updateAction(id, "temperature", inputQuantity);
    } else {
      alert("Please input an integer value");
    }
  };

  return (
    <InputGroup size="sm" className="mb-3" key={id.toString()}>
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">Temperature</InputGroup.Text>
      </InputGroup.Prepend>
      <FormControl
        aria-label="Small"
        aria-describedby="inputGroup-sizing-sm"
        placeholder={Temperature.temperature}
        onChange={(event) => handleTemperatureChange(event)}
      />
    </InputGroup>
  );
};

export default SetTemperature;
