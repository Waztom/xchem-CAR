import React, { useState } from "react";
import axios from "axios";

import InputGroup from "react-bootstrap/InputGroup";
import FormControl from "react-bootstrap/FormControl";

const SetDryingAgent = ({ action, updateAction }) => {
  const dryingagent = action.dryingagent;
  const actiontype = action.actiontype;
  const id = action.id;

  const [DryingAgent, setDryingAgent] = useState({ dryingagent });

  async function patchDryingAgent(value) {
    try {
      const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
        [dryingagent]: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  const handleDryingAgentChange = (e) => {
    const newdryingagent = e.target.value;
    setDryingAgent(newdryingagent);
    patchDryingAgent(newdryingagent);
    updateAction(id, "dryingagent", newdryingagent);
  };

  return (
    <InputGroup size="sm" className="mb-3" key={id.toString()}>
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">
          Drying agent
        </InputGroup.Text>
      </InputGroup.Prepend>
      <FormControl
        aria-label="Small"
        aria-describedby="inputGroup-sizing-sm"
        placeholder={DryingAgent.dryingagent}
        onChange={(event) => handleDryingAgentChange(event)}
      />
    </InputGroup>
  );
};

export default SetDryingAgent;
