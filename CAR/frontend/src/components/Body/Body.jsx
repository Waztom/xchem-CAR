import React, { useState, useEffect } from 'react';
import axios from 'axios';

import Card from 'react-bootstrap/Card';
import ListGroup from 'react-bootstrap/ListGroup';
import Button from 'react-bootstrap/Button';
import InputGroup from 'react-bootstrap/InputGroup';
import FormControl from 'react-bootstrap/FormControl';
import Form from 'react-bootstrap/Form';
import Spinner from 'react-bootstrap/Spinner';

import MethodBody from '../MethodBody/MethodBody';

// Start with main body and then add components
const SetTargetMassInput = ({ targetmass, unit, targetid }) => {
  const [TargetMass, setTargetMass] = useState(targetmass);
  const [Unit, setUnit] = useState(unit);

  async function patchTargetMass(value) {
    try {
      const response = await axios.patch(`api/targets/${targetid}/`, {
        targetmass: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  async function patchTargetUnit(value) {
    try {
      const response = await axios.patch(`api/targets/${targetid}/`, {
        unit: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  const handleTargetMassChange = (e) => {
    const inputTargetMass = e.target.value;

    if (!isNaN(inputTargetMass)) {
      setTargetMass(e.target.value);
      patchTargetMass(Number(e.target.value));
    } else {
      alert('Please input an integer value');
    }
  };

  const handleUnitChange = (e) => {
    setUnit(e.target.value);
    patchTargetUnit(e.target.value);
  };

  return (
    <InputGroup size="sm" className="mb-3">
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">Target Mass</InputGroup.Text>
      </InputGroup.Prepend>
      <FormControl
        aria-label="Small"
        aria-describedby="inputGroup-sizing-sm"
        placeholder={TargetMass}
        onChange={(event) => handleTargetMassChange(event)}
      />
      <Form.Control
        as="select"
        onChange={(event) => handleUnitChange(event)}
        size="sm"
        type="text"
        value={Unit}
      >
        <option>mg</option>
        <option>g</option>
        <option>mmol</option>
      </Form.Control>
    </InputGroup>
  );
};

const TargetCard = ({ target, handleDelete }) => {
  function pressDelete(id) {
    handleDelete(id);
  }

  return (
    <Card text="dark" border="light" style={{ width: '13rem' }}>
      <Card.Header>{target.name}</Card.Header>
      <SetTargetMassInput
        targetmass={target.targetmass}
        unit={target.unit}
        targetid={target.id}
      ></SetTargetMassInput>
      <Card.Body>
        <Card.Img src={target.image} fluid />
      </Card.Body>
      <Button key={target.id} onClick={() => pressDelete(target.id)}>
        Delete Target
      </Button>
    </Card>
  );
};

const Body = ({ ProjectID }) => {
  const [isLoading, setLoading] = useState(true);
  const [Targets, setTargets] = useState([]);

  useEffect(() => {
    async function fetchData() {
      const request = await axios.get(`api/targets?search=${ProjectID}`);
      setTargets(request.data);
      setLoading(false);
    }
    if (ProjectID !== 0) {
      fetchData();
    }
  }, []);

  if (ProjectID !== 0) {
    if (isLoading) {
      return (
        <Spinner animation="border" role="status">
          <span className="sr-only">Loading...</span>
        </Spinner>
      );
    }
  }

  async function deleteData(targetid) {
    try {
      const response = await axios.delete(`api/targets/${targetid}`);
    } catch (error) {
      console.log(error);
    }
  }

  function handleDelete(id) {
    deleteData(id);
    const newList = Targets.filter((item) => item.id !== id);
    setTargets(newList);
  }

  return Targets.map((target) => (
    <ListGroup key={target.id} horizontal>
      <ListGroup.Item key={target.id}>
        <TargetCard target={target} handleDelete={handleDelete} />
      </ListGroup.Item>
      <ListGroup.Item className="targetmethods">
        <MethodBody
          key={target.id}
          targetid={target.id}
          deleteTarget={handleDelete}
        />
      </ListGroup.Item>
    </ListGroup>
  ));
};

export default Body;
