import React, { useState, useEffect } from "react";
import axios from "axios";
import { DragDropContext, Droppable, Draggable } from "react-beautiful-dnd";

import Card from "react-bootstrap/Card";
import ListGroup from "react-bootstrap/ListGroup";
import Button from "react-bootstrap/Button";
import Accordion from "react-bootstrap/Accordion";
import Image from "react-bootstrap/Image";
import CardDeck from "react-bootstrap/CardDeck";
import Spinner from "react-bootstrap/Spinner";
import Container from "react-bootstrap/Container";
import Row from "react-bootstrap/Row";
import Col from "react-bootstrap/Col";
import { Form } from "react-bootstrap";
import InputGroup from "react-bootstrap/InputGroup";
import FormControl from "react-bootstrap/FormControl";

String.prototype.capitalize = function () {
  return this.charAt(0).toUpperCase() + this.slice(1);
};

const SetSolvent = ({ action, updateAction, name }) => {
  const solvname = name + "solvent";

  const solvent = action[solvname];
  const actiontype = action.actiontype;
  const id = action.id;

  const [Solvent, setSolvent] = useState({ solvent });

  async function patchSolvent(value) {
    try {
      const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
        [solvname]: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  const handleSolventChange = (e) => {
    const newsolvent = e.target.value;
    setSolvent(newsolvent);
    patchSolvent(newsolvent);
    updateAction(id, solvname, newsolvent);
  };

  const nameChange = (name) => {
    return name === "rinsingsolvent"
      ? "Rinsing solv"
      : name === "extractionforprecipitatesolvent"
      ? "Extract for precip solv"
      : name === "firstpartitionsolvent"
      ? "1st partition solv"
      : name === "secondpartitionsolvent"
      ? "2nd partition Solv"
      : "Solvent";
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
        placeholder={Solvent.solvent}
        onChange={(event) => handleSolventChange(event)}
      />
    </InputGroup>
  );
};

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

const SetDuration = ({ action, updateAction }) => {
  const duration = action.duration;
  const unit = action.durationunit;
  const actiontype = action.actiontype;
  const id = action.id;

  const [Duration, setDuration] = useState({ duration });
  const [Unit, setUnit] = useState({ unit });

  async function patchDuration(value) {
    try {
      const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
        duration: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  async function patchDurationUnit(value) {
    try {
      const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
        durationunit: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  const handleDurationChange = (e) => {
    const inputQuantity = e.target.value;

    if (!isNaN(inputQuantity)) {
      setDuration(inputQuantity);
      patchDuration(Number(inputQuantity));
      updateAction(id, "duration", inputQuantity);
    } else {
      alert("Please input an integer value");
    }
  };

  const handleUnitChange = (e) => {
    const newunit = e.target.value;
    setUnit(newunit);
    patchDurationUnit(newunit);
    updateAction(id, "durationunit", newunit);
  };

  return (
    <InputGroup size="sm" className="mb-3" key={id.toString()}>
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">Duration</InputGroup.Text>
      </InputGroup.Prepend>
      <FormControl
        aria-label="Small"
        aria-describedby="inputGroup-sizing-sm"
        placeholder={Duration.duration}
        onChange={(event) => handleDurationChange(event)}
      />
      <Form.Control
        as="select"
        onChange={(event) => handleUnitChange(event)}
        size="sm"
        type="text"
        value={Unit.unit}
      >
        <option>seconds</option>
        <option>minutes</option>
        <option>hours</option>
      </Form.Control>
    </InputGroup>
  );
};

const SetStirring = ({ action, updateAction }) => {
  const stirringspeed = action.stirringspeed;
  const actiontype = action.actiontype;
  const id = action.id;

  const [StirringSpeed, setStirringSpeed] = useState({ stirringspeed });

  async function patchStirringSpeed(value) {
    try {
      const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
        stirringspeed: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  const handleStirringSpeedChange = (e) => {
    const newspeed = e.target.value;
    setStirringSpeed(newspeed);
    patchStirringSpeed(nwspeed);
    updateAction(id, "stirringspeed", newspeed);
  };

  return (
    <InputGroup size="sm" className="mb-3" key={id.toString()}>
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">
          Stirring speed
        </InputGroup.Text>
      </InputGroup.Prepend>
      <Form.Control
        as="select"
        onChange={(event) => handleStirringSpeedChange(event)}
        size="sm"
        type="text"
        value={StirringSpeed.stirringspeed}
      >
        <option>gentle</option>
        <option>normal</option>
        <option>vigourous</option>
      </Form.Control>
    </InputGroup>
  );
};

const SetTemperature = ({ action, updateAction }) => {
  const temperature = action.temperature;
  const actiontype = action.actiontype;
  const id = action.id;

  const [Temperature, SetTemperature] = useState({ temperature });

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
      SetTemperature(inputQuantity);
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

const SetNumberRepetitions = ({ action, updateAction }) => {
  const numberrepetitions = action.numberofrepetitions;
  const actiontype = action.actiontype;
  const id = action.id;

  const [NumberRepetitions, SetNumberRepetitions] = useState({
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
      SetNumberRepetitions(inputQuantity);
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

const SetpH = ({ action, updateAction }) => {
  const ph = action.pH;
  const actiontype = action.actiontype;
  const id = action.id;

  const [pH, setpH] = useState({ ph });

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
        placeholder={pH.ph}
        onChange={(event) => handlepHChange(event)}
      />
    </InputGroup>
  );
};

const SetMaterial = ({ action, updateAction }) => {
  const material = action.material;
  const actiontype = action.actiontype;
  const id = action.id;

  const [Material, setMaterial] = useState({ material });

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
        placeholder={Material.material}
        onChange={(event) => handleMaterialChange(event)}
      />
    </InputGroup>
  );
};

const SetGas = ({ action, updateAction }) => {
  const gas = action.gas;
  const actiontype = action.actiontype;
  const id = action.id;

  const [Gas, setGas] = useState({ gas });

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
        placeholder={Gas.gas}
        onChange={(event) => handleGasChange(event)}
      />
    </InputGroup>
  );
};

const SetLayer = ({ action, updateAction }) => {
  const layer = action.layer.capitalize();
  const actiontype = action.actiontype;
  const id = action.id;

  const [Layer, SetLayer] = useState({ layer });

  async function patchLayer(value) {
    try {
      const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
        layer: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  const handleLayerChange = (e) => {
    const newlayer = e.target.value.toLowerCase();
    SetLayer(e.target.value);
    patchLayer(newlayer);
    updateAction(id, "layer", newlayer);
  };

  return (
    <InputGroup size="sm" className="mb-3" key={id.toString()}>
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">Layer</InputGroup.Text>
      </InputGroup.Prepend>
      <Form.Control
        as="select"
        onChange={(event) => handleLayerChange(event)}
        size="sm"
        type="text"
        value={Layer.layer}
      >
        <option>Organic</option>
        <option>Aqueous</option>
      </Form.Control>
    </InputGroup>
  );
};

const SetAtmosphere = ({ action, updateAction }) => {
  const atmosphere = action.atmosphere.capitalize();
  const actiontype = action.actiontype;
  const id = action.id;

  const [Atmosphere, SetAtmosphere] = useState({ atmosphere });

  async function patchAtmosphere(value) {
    try {
      const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
        atmosphere: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  const handleAtmosphereChange = (e) => {
    const newatmosphere = e.target.value.toLowerCase();
    SetAtmosphere(e.target.value);
    patchAtmosphere(newatmosphere);
    updateAction(id, "atmosphere", newatmosphere);
  };

  return (
    <InputGroup size="sm" className="mb-3" key={id.toString()}>
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">Atmosphere</InputGroup.Text>
      </InputGroup.Prepend>
      <Form.Control
        as="select"
        onChange={(event) => handleAtmosphereChange(event)}
        size="sm"
        type="text"
        value={Atmosphere.atmosphere}
      >
        <option>Nitrogen</option>
        <option>Air</option>
      </Form.Control>
    </InputGroup>
  );
};

const SetDropwise = ({ action, updateAction }) => {
  const dropwise = action.dropwise;
  const actiontype = action.actiontype;
  const id = action.id;

  const [DropWise, SetDropwise] = useState({ dropwise });

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
    SetDropwise(newdropwise);
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
        value={DropWise.dropwise}
      >
        <option>True</option>
        <option>False</option>
      </Form.Control>
    </InputGroup>
  );
};

const SetDeanStark = ({ action, updateAction }) => {
  const deanstark = action.deanstarkapparatus;
  const actiontype = action.actiontype;
  const id = action.id;

  const [DeanStark, SetDeanStark] = useState({ deanstark });

  async function patchDeanStark(value) {
    try {
      const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
        deanstarkapparatus: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  const handleDeanStarkChange = (e) => {
    const newdeanstark = e.target.value;
    SetDeanStark(newdeanstark);
    patchDeanStark(newdeanstark);
    updateAction(id, "deanstarkapparatus", newdeanstark);
  };

  return (
    <InputGroup size="sm" className="mb-3" key={id.toString()}>
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">Dean Stark</InputGroup.Text>
      </InputGroup.Prepend>
      <Form.Control
        as="select"
        onChange={(event) => handleDeanStarkChange(event)}
        size="sm"
        type="text"
        value={DeanStark.deanstark}
      >
        <option>True</option>
        <option>False</option>
      </Form.Control>
    </InputGroup>
  );
};

const SetQuantityInput = ({ action, updateAction, name }) => {
  const quantname = name + "quantity";
  const unitname = name + "quantityunit";

  const quantity = action[quantname];
  const unit = action[unitname];
  const actiontype = action.actiontype;
  const id = action.id;

  const [Quantity, setQuantity] = useState({ quantity });
  const [Unit, setUnit] = useState({ unit });

  async function patchQuantity(value) {
    try {
      const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
        [quantname]: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  async function patchQuantityUnit(value) {
    try {
      const response = await axios.patch(`api/IBM${actiontype}actions/${id}/`, {
        [unitname]: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  const handleQuantityChange = (e) => {
    const inputQuantity = e.target.value;

    if (!isNaN(inputQuantity)) {
      setQuantity(e.target.value);
      patchQuantity(Number(e.target.value));
      updateAction(id, quantname, Number(inputQuantity));
    } else {
      alert("Please input an integer value");
    }
  };

  const handleUnitChange = (e) => {
    const inputUnit = e.target.value;

    setUnit(inputUnit);
    patchQuantityUnit(inputUnit);
    updateAction(id, unitname, inputUnit);
  };

  return (
    <InputGroup size="sm" className="mb-3" key={id.toString()}>
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">Quantity</InputGroup.Text>
      </InputGroup.Prepend>
      <FormControl
        aria-label="Small"
        aria-describedby="inputGroup-sizing-sm"
        placeholder={Quantity.quantity}
        onChange={(event) => handleQuantityChange(event)}
      />
      <Form.Control
        as="select"
        onChange={(event) => handleUnitChange(event)}
        size="sm"
        type="text"
        value={Unit.unit}
      >
        <option>moleq</option>
        <option>ml</option>
        <option>mmol</option>
      </Form.Control>
    </InputGroup>
  );
};

const IBMAddAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container key={actionno}>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <Image
            width={150}
            height={150}
            src={action.materialimage}
            alt={action.material}
            fluid
          />
          <SetQuantityInput
            action={action}
            updateAction={updateAction}
            name={"material"}
          ></SetQuantityInput>
          <SetDropwise
            action={action}
            updateAction={updateAction}
          ></SetDropwise>
          <SetAtmosphere
            action={action}
            updateAction={updateAction}
          ></SetAtmosphere>
        </Col>
      </Row>
    </Container>
  );
};

const IBMCollectLayerAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetLayer action={action} updateAction={updateAction}></SetLayer>
        </Col>
      </Row>
    </Container>
  );
};

const IBMConcentrateAction = ({ action, actionno }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <h5>
      {actionno}. {actiontype}
    </h5>
  );
};

const IBMDegasAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetGas action={action} updateAction={updateAction}></SetGas>
          <SetDuration
            action={action}
            updateAction={updateAction}
          ></SetDuration>
        </Col>
      </Row>
    </Container>
  );
};

const IBMDrySolidAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetTemperature
            action={action}
            updateAction={updateAction}
          ></SetTemperature>
          <SetDuration
            action={action}
            updateAction={updateAction}
          ></SetDuration>
          <SetAtmosphere
            action={action}
            updateAction={updateAction}
          ></SetAtmosphere>
        </Col>
      </Row>
    </Container>
  );
};

const IBMDrySolutionAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetDryingAgent
            action={action}
            updateAction={updateAction}
          ></SetDryingAgent>
        </Col>
      </Row>
    </Container>
  );
};

const IBMExtractAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetSolvent action={action} updateAction={updateAction}></SetSolvent>
          <SetQuantityInput
            action={action}
            updateAction={updateAction}
            name={"solvent"}
          ></SetQuantityInput>
          <SetNumberRepetitions
            action={action}
            updateAction={updateAction}
          ></SetNumberRepetitions>
        </Col>
      </Row>
    </Container>
  );
};

const IBMFilterAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetSolvent
            action={action}
            updateAction={updateAction}
            name={"rinsing"}
          ></SetSolvent>
          <SetQuantityInput
            action={action}
            updateAction={updateAction}
            name={"rinsingsolvent"}
          ></SetQuantityInput>
          <SetSolvent
            action={action}
            updateAction={updateAction}
            name={"extractionforprecipitate"}
          ></SetSolvent>
          <SetQuantityInput
            action={action}
            updateAction={updateAction}
            name={"extractionforprecipitatesolvent"}
          ></SetQuantityInput>
        </Col>
      </Row>
    </Container>
  );
};

const IBMMakeSolutionAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <Image
            width={150}
            height={150}
            src={action.soluteimage}
            alt={action.solute}
            fluid
          />
          <SetQuantityInput
            action={action}
            updateAction={updateAction}
            name={"solute"}
          ></SetQuantityInput>
          <Image
            width={150}
            height={150}
            src={action.solventimage}
            alt={action.solvent}
            fluid
          />
          <SetQuantityInput
            action={action}
            updateAction={updateAction}
            name={"solvent"}
          ></SetQuantityInput>
        </Col>
      </Row>
    </Container>
  );
};

const IBMPartitionAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetSolvent
            action={action}
            updateAction={updateAction}
            name={"firstparition"}
          ></SetSolvent>
          <SetQuantityInput
            action={action}
            updateAction={updateAction}
            name={"firstpartitionsolvent"}
          ></SetQuantityInput>
          <SetSolvent
            action={action}
            updateAction={updateAction}
            name={"secondpartition"}
          ></SetSolvent>
          <SetQuantityInput
            action={action}
            updateAction={updateAction}
            name={"secondpartitionsolvent"}
          ></SetQuantityInput>
        </Col>
      </Row>
    </Container>
  );
};

const IBMpHAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container key={actionno}>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetMaterial
            action={action}
            updateAction={updateAction}
          ></SetMaterial>
          <SetQuantityInput
            action={action}
            updateAction={updateAction}
            name={"material"}
          ></SetQuantityInput>
          <SetpH action={action} updateAction={updateAction}></SetpH>
          <SetDropwise
            action={action}
            updateAction={updateAction}
          ></SetDropwise>
          <SetTemperature
            action={action}
            updateAction={updateAction}
          ></SetTemperature>
        </Col>
      </Row>
    </Container>
  );
};

const IBMPhaseSeparationAction = ({ action, actionno }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <h5>
      {actionno}. {actiontype}
    </h5>
  );
};

const IBMQuenchAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container key={actionno}>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetMaterial
            action={action}
            updateAction={updateAction}
          ></SetMaterial>
          <SetQuantityInput
            action={action}
            updateAction={updateAction}
            name={"material"}
          ></SetQuantityInput>
          <SetDropwise
            action={action}
            updateAction={updateAction}
          ></SetDropwise>
          <SetTemperature
            action={action}
            updateAction={updateAction}
          ></SetTemperature>
        </Col>
      </Row>
    </Container>
  );
};

const IBMRefluxAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container key={actionno}>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetDuration
            action={action}
            updateAction={updateAction}
          ></SetDuration>

          <SetStirring
            action={action}
            updateAction={updateAction}
          ></SetStirring>
          <SetDeanStark
            action={action}
            updateAction={updateAction}
          ></SetDeanStark>

          <SetAtmosphere
            action={action}
            updateAction={updateAction}
          ></SetAtmosphere>
        </Col>
      </Row>
    </Container>
  );
};

const IBMSetTemperatureAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container key={actionno}>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetTemperature
            action={action}
            updateAction={updateAction}
          ></SetTemperature>
        </Col>
      </Row>
    </Container>
  );
};

const IBMStirAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container key={actionno}>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetDuration
            action={action}
            updateAction={updateAction}
          ></SetDuration>
          <SetStirring
            action={action}
            updateAction={updateAction}
          ></SetStirring>
          <SetAtmosphere
            action={action}
            updateAction={updateAction}
          ></SetAtmosphere>
          <SetTemperature
            action={action}
            updateAction={updateAction}
          ></SetTemperature>
        </Col>
      </Row>
    </Container>
  );
};

const IBMStoreAction = ({ action, actionno }) => {
  const actiontype = action.actiontype.capitalize();
  const material = action.material.capitalize();

  return (
    <Container key={actionno}>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <Image
            width={150}
            height={150}
            src={action.materialimage}
            alt={action.material}
            fluid
          />
        </Col>
      </Row>
    </Container>
  );
};

const IBMWaitAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container key={actionno}>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetDuration
            action={action}
            updateAction={updateAction}
          ></SetDuration>
          <SetTemperature
            action={action}
            updateAction={updateAction}
          ></SetTemperature>
        </Col>
      </Row>
    </Container>
  );
};

const IBMWashAction = ({ action, actionno, updateAction }) => {
  const actiontype = action.actiontype.capitalize();

  return (
    <Container key={actionno}>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <SetMaterial
            action={action}
            updateAction={updateAction}
          ></SetMaterial>
          <SetQuantityInput
            action={action}
            updateAction={updateAction}
            name={"material"}
          ></SetQuantityInput>
          <SetNumberRepetitions
            action={action}
            updateAction={updateAction}
          ></SetNumberRepetitions>
        </Col>
      </Row>
    </Container>
  );
};

const ActionsList = ({ reactionid }) => {
  const [isLoading, setLoading] = useState(true);
  const [Actions, setActions] = useState([]);

  function compareActionno(a, b) {
    return a.actionno - b.actionno;
  }

  useEffect(() => {
    async function fetchData() {
      const URLs = [
        `api/IBMaddactions?search=${reactionid}`,
        `api/IBMcollect-layeractions?search=${reactionid}`,
        `api/IBMconcentrateactions?search=${reactionid}`,
        `api/IBMdegasactions?search=${reactionid}`,
        `api/IBMdry-solidactions?search=${reactionid}`,
        `api/IBMdry-solutionactions?search=${reactionid}`,
        `api/IBMextractactions?search=${reactionid}`,
        `api/IBMfilteractions?search=${reactionid}`,
        `api/IBMmake-solutionactions?search=${reactionid}`,
        `api/IBMpartitionactions?search=${reactionid}`,
        `api/IBMphactions?search=${reactionid}`,
        `api/IBMphase-separationactions?search=${reactionid}`,
        `api/IBMquenchactions?search=${reactionid}`,
        `api/IBMrefluxactions?search=${reactionid}`,
        `api/IBMset-temperatureactions?search=${reactionid}`,
        `api/IBMstiractions?search=${reactionid}`,
        `api/IBMstoreactions?search=${reactionid}`,
        `api/IBMwaitactions?search=${reactionid}`,
        `api/IBMwashactions?search=${reactionid}`,
      ];
      const requests = URLs.map((URL) => axios.get(URL).catch((err) => null));

      try {
        const [
          addresp,
          collectresp,
          concresp,
          degasresp,
          drysolidresp,
          drysolnresp,
          extractresp,
          filterresp,
          makesolnresp,
          partitionresp,
          phresp,
          phasesepresp,
          quenchresp,
          refluxresp,
          settempresp,
          stirresp,
          storeresp,
          waitresp,
          washresp,
        ] = await axios.all(requests);
        const actions = addresp.data.concat(
          collectresp.data,
          concresp.data,
          degasresp.data,
          drysolidresp.data,
          drysolnresp.data,
          extractresp.data,
          filterresp.data,
          makesolnresp.data,
          partitionresp.data,
          phresp.data,
          phasesepresp.data,
          quenchresp.data,
          refluxresp.data,
          settempresp.data,
          stirresp.data,
          storeresp.data,
          waitresp.data,
          washresp.data
        );
        const sorted = actions.sort(compareActionno);
        setActions(sorted);
        setLoading(false);
      } catch (err) {
        console.log(err);
      }
    }
    fetchData();
  }, []);

  if (isLoading) {
    return (
      <Spinner animation="border" role="status">
        <span className="sr-only">Loading...</span>
      </Spinner>
    );
  }

  async function patchActionNo(actiontype, actionid, value) {
    try {
      const response = await axios.patch(
        `api/IBM${actiontype}actions/${actionid}/`,
        {
          actionno: value,
        }
      );
    } catch (error) {
      console.log(error);
    }
  }

  function handleOnDragEnd(result) {
    const items = Array.from(Actions);
    const [reordereditem] = items.splice(result.source.index, 1);
    items.splice(result.destination.index, 0, reordereditem);

    setActions(items);

    items.map((item, index) => {
      const actiontype = item.actiontype;
      const newactionno = index + 1;
      patchActionNo(actiontype, item.id, newactionno);
    });
  }

  function getComponent(action, actionno, updateAction) {
    const components = {
      add: (
        <IBMAddAction
          action={action}
          actionno={actionno}
          updateAction={updateAction}
        ></IBMAddAction>
      ),
      "collect-layer": (
        <IBMCollectLayerAction
          action={action}
          actionno={actionno}
          updateAction={updateAction}
        ></IBMCollectLayerAction>
      ),
      concentrate: (
        <IBMConcentrateAction
          action={action}
          actionno={actionno}
        ></IBMConcentrateAction>
      ),
      degas: (
        <IBMDegasAction
          action={action}
          actionno={actionno}
          updateAction={updateAction}
        ></IBMDegasAction>
      ),
      "dry-solid": (
        <IBMDrySolidAction
          action={action}
          actionno={actionno}
          updateAction={updateAction}
        ></IBMDrySolidAction>
      ),
      "dry-solution": (
        <IBMDrySolutionAction
          action={action}
          actionno={actionno}
          updateAction={updateAction}
        ></IBMDrySolutionAction>
      ),
      extract: (
        <IBMExtractAction
          action={action}
          actionno={actionno}
          updateAction={updateAction}
        ></IBMExtractAction>
      ),
      filter: (
        <IBMFilterAction
          action={action}
          actionno={actionno}
          updateAction={updateAction}
        ></IBMFilterAction>
      ),
      "make-solution": (
        <IBMMakeSolutionAction
          action={action}
          actionno={actionno}
          updateAction={updateAction}
        ></IBMMakeSolutionAction>
      ),
      partition: (
        <IBMPartitionAction
          action={action}
          actionno={actionno}
          updateAction={updateAction}
        ></IBMPartitionAction>
      ),
      ph: (
        <IBMpHAction
          action={action}
          actionno={actionno}
          updateAction={updateAction}
        ></IBMpHAction>
      ),
      "phase-separation": (
        <IBMPhaseSeparationAction
          action={action}
          actionno={actionno}
          updateAction={updateAction}
        ></IBMPhaseSeparationAction>
      ),
      quench: (
        <IBMQuenchAction
          action={action}
          actionno={actionno}
          updateAction={updateAction}
        ></IBMQuenchAction>
      ),
      reflux: (
        <IBMRefluxAction
          action={action}
          actionno={actionno}
          updateAction={updateAction}
        ></IBMRefluxAction>
      ),
      "set-temperature": (
        <IBMSetTemperatureAction
          action={action}
          actionno={actionno}
          updateAction={updateAction}
        ></IBMSetTemperatureAction>
      ),
      stir: (
        <IBMStirAction
          action={action}
          actionno={actionno}
          updateAction={updateAction}
        ></IBMStirAction>
      ),
      store: (
        <IBMStoreAction
          action={action}
          actionno={actionno}
          updateAction={updateAction}
        ></IBMStoreAction>
      ),
      wait: (
        <IBMWaitAction
          action={action}
          actionno={actionno}
          updateAction={updateAction}
        ></IBMWaitAction>
      ),
      wash: (
        <IBMWashAction
          action={action}
          actionno={actionno}
          updateAction={updateAction}
        ></IBMWashAction>
      ),
    };

    return components[action.actiontype];
  }

  function updateAction(actionid, changekey, changevalue) {
    // let newArr = [...Actions];
    var foundIndex = Actions.findIndex((x) => x.id == actionid);
    Actions[foundIndex][changekey] = changevalue;
    // setActions(newArr);
    console.log(Actions[foundIndex]);
  }

  return (
    <DragDropContext onDragEnd={handleOnDragEnd}>
      <Droppable droppableId="characters">
        {(provided) => (
          <ListGroup
            className="characters"
            {...provided.droppableProps}
            ref={provided.innerRef}
          >
            {Actions.map((action, index) => {
              const actionid =
                action.actiontype +
                "-" +
                action.id.toString() +
                "-" +
                index.toString();
              const actionno = index + 1;
              return (
                <Draggable key={actionid} draggableId={actionid} index={index}>
                  {(provided) => (
                    <ListGroup.Item
                      ref={provided.innerRef}
                      {...provided.draggableProps}
                      {...provided.dragHandleProps}
                    >
                      {getComponent(action, actionno, updateAction)}
                    </ListGroup.Item>
                  )}
                </Draggable>
              );
            })}
            {provided.placeholder}
          </ListGroup>
        )}
      </Droppable>
    </DragDropContext>
  );
};

const ProductImage = ({ reactionid }) => {
  // Use hooks instead of classes
  const [isLoading, setLoading] = useState(true);
  const [Product, setProduct] = useState([]);

  useEffect(() => {
    async function fetchData() {
      try {
        const request = await axios.get(`api/products?search=${reactionid}`);
        setProduct(request.data);
        setLoading(false);
      } catch (err) {
        console.log(err);
      }
    }
    fetchData();
  }, []);

  if (isLoading) {
    return (
      <Spinner animation="border" role="status">
        <span className="sr-only">Loading...</span>
      </Spinner>
    );
  }

  return Product.map((product) => <Image src={product.image} fluid />);
};

const ReactionAccordian = ({ methodid }) => {
  // Use hooks instead of classes
  const [isLoading, setLoading] = useState(true);
  const [Reactions, setReactions] = useState([]);

  useEffect(() => {
    async function fetchData() {
      try {
        const request = await axios.get(`api/reactions?search=${methodid}`);
        setReactions(request.data);
        setLoading(false);
      } catch (err) {
        console.log(err);
      }
    }
    fetchData();
  }, []);

  if (isLoading) {
    return (
      <Spinner animation="border" role="status">
        <span className="sr-only">Loading...</span>
      </Spinner>
    );
  }

  return (
    <Accordion>
      {Reactions.map((reaction) => (
        <Card key={reaction.id}>
          <Card.Header>
            <Accordion.Toggle as={Button} variant="link" eventKey={reaction.id}>
              <ProductImage
                key={reaction.id}
                reactionid={reaction.id}
              ></ProductImage>
              {reaction.reactionclass}
            </Accordion.Toggle>
          </Card.Header>
          <Accordion.Collapse eventKey={reaction.id}>
            <ActionsList
              key={reaction.id}
              reactionid={reaction.id}
            ></ActionsList>
          </Accordion.Collapse>
        </Card>
      ))}
    </Accordion>
  );
};

const MethodCard = ({ method }) => {
  return (
    <CardDeck>
      <Card border="light" style={{ width: "30rem" }} key={method.id}>
        <Card.Body>
          <Card.Title>Synthetic steps</Card.Title>
          <ReactionAccordian
            key={method.id}
            methodid={method.id}
          ></ReactionAccordian>
        </Card.Body>
      </Card>
    </CardDeck>
  );
};

const MethodBody = ({ targetid, deleteTarget }) => {
  // Use hooks instead of classes
  const [isLoading, setLoading] = useState(true);
  const [Methods, setMethods] = useState([]);
  const [NoMethods, setNoMethods] = useState([]);

  useEffect(() => {
    async function fetchData() {
      try {
        const request = await axios.get(`api/methods?search=${targetid}`);
        setMethods(request.data);
        const nomethods = Object.keys(request.data).length;
        setNoMethods(nomethods);
        setLoading(false);
      } catch (err) {}
    }
    fetchData();
  }, []);

  if (isLoading) {
    return (
      <Spinner animation="border" role="status">
        <span className="sr-only">Loading...</span>
      </Spinner>
    );
  }

  async function deleteData(methodid) {
    try {
      const response = await axios.delete(`api/methods/${methodid}`);
    } catch (error) {
      console.log(error);
    }
  }

  function handleDelete(id) {
    deleteData(id);
    const newList = Methods.filter((item) => item.id !== id);
    setMethods(newList);
    const nomethods = Object.keys(Methods).length;
    if (nomethods == 1) {
      deleteTarget(targetid);
    }
  }

  if (NoMethods == 0) {
    return (
      <Card bg="warning">
        <Card.Body>
          Backend still busy processing else no methods were found!
        </Card.Body>
      </Card>
    );
  } else {
    return (
      <ListGroup horizontal>
        {Methods.map((method) => (
          <ListGroup.Item key={method.id}>
            <MethodCard method={method} key={method.id} />
            <Button key={method.id} onClick={() => handleDelete(method.id)}>
              Delete Method
            </Button>
          </ListGroup.Item>
        ))}
      </ListGroup>
    );
  }
};

export default MethodBody;
