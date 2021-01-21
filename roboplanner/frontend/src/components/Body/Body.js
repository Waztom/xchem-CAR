import React, { useState, useEffect } from "react";
import axios from "axios";

import Card from "react-bootstrap/Card";
import ListGroup from "react-bootstrap/ListGroup";
import Button from "react-bootstrap/Button";

import MethodBody from "../MethodBody/MethodBody";

// Start with main body and then add components
const TargetCard = ({ name, image }) => {
  return (
    <Card text="dark" border="light" style={{ width: "18rem" }}>
      <Card.Header>{name}</Card.Header>
      <Card.Body>
        <Card.Img variant="bottom" src={image} />
      </Card.Body>
    </Card>
  );
};

const Body = ({ ProjectID }) => {
  // Use hooks instead of classes
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
      return <div className="App">Loading...</div>;
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
        <TargetCard name={target.name} image={target.image} />
        <Button key={target.id} onClick={() => handleDelete(target.id)}>
          Delete Target
        </Button>
      </ListGroup.Item>
      <ListGroup.Item>
        <MethodBody
          key={target.id}
          targetid={target.id}
          onChange={handleDelete}
        />
      </ListGroup.Item>
    </ListGroup>
  ));
};

export default Body;
