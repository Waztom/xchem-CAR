import React, { useState, useEffect } from "react";
import axios from "axios";

import Card from "react-bootstrap/Card";
import ListGroup from "react-bootstrap/ListGroup";

import ReactionBody from "../ReactionBody/ReactionBody";

// Start with main body and then add components
const TargetCard = ({ name, image }) => {
  return (
    <Card style={{ width: "18rem" }}>
      <Card.Img variant="top" src={image} />
      <Card.Body>
        <Card.Title>{name}</Card.Title>
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
    if (ProjectID !== undefined || ProjectID !== 0) {
      fetchData();
    }
  }, []);

  if (ProjectID !== undefined || ProjectID !== 0) {
    if (isLoading) {
      return <div className="App">Loading...</div>;
    }
  }

  return Targets.map((target) => (
    <ListGroup horizontal key={target.name}>
      <ListGroup.Item key={target.name}>
        <TargetCard key={target.name} name={target.name} image={target.image} />
      </ListGroup.Item>
      {/* <ListGroup.Item>
        <ReactionBody targetid={target.id} />
      </ListGroup.Item> */}
    </ListGroup>

    // <React.Fragment>
    //   <Card>
    //     <Container fluid>
    //       <Row>
    //         <Col sm={3}>
    //           <TargetCard
    //             key={target.uniqueId}
    //             name={target.name}
    //             image={target.image}
    //           />
    //         </Col>
    //         <Col sm={8}>
    //           <ReactionBody targetid={target.id} />
    //         </Col>
    //       </Row>
    //     </Container>
    //   </Card>
    // </React.Fragment>
  ));
};

export default Body;
