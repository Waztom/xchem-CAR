import React, { useState } from 'react';

import InputGroup from 'react-bootstrap/InputGroup';
import FormControl from 'react-bootstrap/FormControl';
import { patchChange } from '../Utils';

const SetSolvent = ({ action, updateAction, name }) => {
  const checkName = (name) => {
    return name === 'rinsing'
      ? 'rinsingsolvent'
      : name === 'extractionforprecipitate'
      ? 'extractionforprecipitatesolvent'
      : name === 'firstpartition'
      ? 'firstpartitionsolvent'
      : name === 'secondpartition'
      ? 'secondpartitionsolvent'
      : 'solvent';
  };

  const solvname = checkName(name);
  const solvent = action[solvname];
  const actiontype = action.actiontype;
  const id = action.id;

  const [Solvent, setSolvent] = useState(solvent);

  const handleSolventChange = (e) => {
    console.log(action);
    console.log(solvname);

    const newsolvent = e.target.value;
    setSolvent(newsolvent);
    patchChange(actiontype, id, solvname, newsolvent);
    updateAction(id, solvname, newsolvent);
  };

  const nameChange = (name) => {
    return name === 'rinsingsolvent'
      ? 'Rinsing solv'
      : name === 'extractionforprecipitatesolvent'
      ? 'Extract for precip solv'
      : name === 'firstpartitionsolvent'
      ? '1st partition solv'
      : name === 'secondpartitionsolvent'
      ? '2nd partition Solv'
      : 'Solvent';
  };

  return (
    <InputGroup size="sm" className="mb-3" key={id.toString()}>
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">
          {nameChange(solvname)}
        </InputGroup.Text>
      </InputGroup.Prepend>
      <FormControl
        aria-label="Small"
        aria-describedby="inputGroup-sizing-sm"
        placeholder={Solvent}
        onChange={(event) => handleSolventChange(event)}
      />
    </InputGroup>
  );
};

export default SetSolvent;
