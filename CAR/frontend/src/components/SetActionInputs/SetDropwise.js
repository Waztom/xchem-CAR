import React, { useState } from "react";
import axios from "axios";

import { Form } from "react-bootstrap";
import InputGroup from "react-bootstrap/InputGroup";

const SetDropwise = ({ action, updateAction }) => {
  const dropwise = action.dropwise;
  const actiontype = action.actiontype;
  const id = action.id;

  const [DropWise, setDropwise] = useState(dropwise);

  async function patchDropWise(value) {
    try {
      const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
        dropwise: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  const handleDropWiseChange = (e) => {
    const newdropwise = e.target.value;
    setDropwise(newdropwise);
    patchDropWise(newdropwise);
    updateAction(id, "dropwise", newdropwise);
  };

  return (
    <InputGroup size="sm" className="mb-3" key={id.toString()}>
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">Dropwise</InputGroup.Text>
      </InputGroup.Prepend>
      <Form.Control
        as="select"
        onChange={(event) => handleDropWiseChange(event)}
        size="sm"
        type="text"
        value={DropWise}
      >
        <option>True</option>
        <option>False</option>
      </Form.Control>
    </InputGroup>
  );
};

export default SetDropwise;
