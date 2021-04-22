import React, { useState } from "react";
import axios from "axios";

import InputGroup from "react-bootstrap/InputGroup";
import FormControl from "react-bootstrap/FormControl";

const SetMaterial = ({ action, updateAction }) => {
  const material = action.material;
  const actiontype = action.actiontype;
  const id = action.id;

  const [Material, setMaterial] = useState(material);

  async function patchMaterial(value) {
    try {
      const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
        material: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  const handleMaterialChange = (e) => {
    const newmaterial = e.target.value;
    setMaterial(newmaterial);
    patchMaterial(newmaterial);
    updateAction(id, "material", newmaterial);
  };

  return (
    <InputGroup size="sm" className="mb-3" key={id.toString()}>
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">Material</InputGroup.Text>
      </InputGroup.Prepend>
      <FormControl
        aria-label="Small"
        aria-describedby="inputGroup-sizing-sm"
        placeholder={Material}
        onChange={(event) => handleMaterialChange(event)}
      />
    </InputGroup>
  );
};

export default SetMaterial;
