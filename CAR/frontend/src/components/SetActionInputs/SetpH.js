import React, { useState } from "react";
import axios from "axios";

import InputGroup from "react-bootstrap/InputGroup";
import FormControl from "react-bootstrap/FormControl";

const SetpH = ({ action, updateAction }) => {
  const ph = action.pH;
  const actiontype = action.actiontype;
  const id = action.id;

  const [pH, setpH] = useState(ph);

  async function patchpH(value) {
    try {
      const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
        pH: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  const handlepHChange = (e) => {
    const inputQuantity = e.target.value;

    if (!isNaN(inputQuantity)) {
      setpH(inputQuantity);
      patchpH(Number(inputQuantity));
      updateAction(id, "pH", inputQuantity);
    } else {
      alert("Please input an integer value");
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
      />
    </InputGroup>
  );
};

export default SetpH;
