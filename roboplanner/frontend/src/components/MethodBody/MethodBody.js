import React, { useState, useEffect } from "react";
import axios from "axios";
import { DragDropContext, Droppable, Draggable } from "react-beautiful-dnd";

import Card from "react-bootstrap/Card";
import ListGroup from "react-bootstrap/ListGroup";
import Container from "react-bootstrap/Container";
import Col from "react-bootstrap/Col";
import Row from "react-bootstrap/Row";
import InputGroup from "react-bootstrap/InputGroup";
import FormControl from "react-bootstrap/FormControl";
import Form from "react-bootstrap/Form";
import Button from "react-bootstrap/Button";
import Accordion from "react-bootstrap/Accordion";
import Image from "react-bootstrap/Image";
import CardDeck from "react-bootstrap/CardDeck";
import { List } from "react-bootstrap-icons";

//Start with most nested children components
// Process cards - these cards describe the overall
// process: order of reactant/reagent addition, work-up,
// QC etc. Will define these actions and breakpoints as
// we go

const SetMoleqInput = ({ moleq }) => {
  const [Moleq, setMoleq] = useState({ moleq });

  const handleMoleqChange = (e) => {
    const inputMoleq = e.target.value;

    if (!isNaN(inputMoleq)) {
      setMoleq(e.target.value);
    } else {
      alert("Please input an integer value");
    }

    // Add code to update API as well
  };

  return (
    <InputGroup
      size="sm"
      className="mb-3"
      value={Moleq.moleq}
      onChange={handleMoleqChange}
    >
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">Mol eq.</InputGroup.Text>
      </InputGroup.Prepend>
      <FormControl
        aria-label="Small"
        aria-describedby="inputGroup-sizing-sm"
        placeholder={Moleq.moleq}
      />
      {console.log(Moleq)}
    </InputGroup>
  );
};

const SetAddOrderInput = ({ addorder }) => {
  const [AddOrder, setAddOrder] = useState({ addorder });

  const handleAddOrderChange = (e) => {
    const inputOrder = e.target.value;

    if (!isNaN(inputOrder)) {
      setAddOrder(e.target.value);
    } else {
      alert("Please input an integer value");
    }

    // Add code to update API as well
  };

  return (
    <InputGroup
      size="sm"
      className="mb-3"
      value={AddOrder.addorder}
      onChange={handleAddOrderChange}
    >
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">Add order</InputGroup.Text>
      </InputGroup.Prepend>
      <FormControl
        aria-label="Small"
        aria-describedby="inputGroup-sizing-sm"
        placeholder={AddOrder.addorder}
      />
      {console.log(AddOrder)}
    </InputGroup>
  );
};

const SetConcentrationInput = ({ concentration }) => {
  const [Concentration, setConcentration] = useState({ concentration });

  const handleConcentrationChange = (e) => {
    const inputConcentration = e.target.value;

    if (!isNaN(inputConcentration)) {
      setConcentration(e.target.value);
    } else {
      alert("Please input an integer value");
    }

    // Add code to update API as well
  };

  return (
    <InputGroup
      size="sm"
      className="mb-3"
      value={Concentration.concentration}
      onChange={handleConcentrationChange}
    >
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">Conc (M)</InputGroup.Text>
      </InputGroup.Prepend>
      <FormControl
        aria-label="Small"
        aria-describedby="inputGroup-sizing-sm"
        placeholder={Concentration.concentration}
      />
      {console.log(Concentration)}
    </InputGroup>
  );
};

const SetSolventInput = ({ solvent }) => {
  const [Solvent, setSolvent] = useState({ solvent });

  const handleSolventChange = (e) => {
    setSolvent(e.target.value);
    // Add code to update API as well
  };

  return (
    <Form>
      <Form.Group controlId="exampleForm.ControlSelect1">
        <Form.Control
          as="select"
          onChange={handleSolventChange}
          size="sm"
          type="text"
        >
          <option>{Solvent.solvent}</option>
          <option>DMA</option>
          <option>DCM</option>
        </Form.Control>
      </Form.Group>
      {console.log(Solvent)}
    </Form>
  );
};

// Reactant/Reactant cards
const ReactantCard = ({ image, moleq, addorder, solvent, concentration }) => {
  // Set state

  // Render to return
  return (
    <Card style={{ width: "12rem" }}>
      {/* Reactant image here */}
      <Card.Img variant="top" src={image} alt="Reactant Image" />
      <Card.Body></Card.Body>
      <Container>
        <Row>
          <Col>
            <SetMoleqInput moleq={moleq} />
          </Col>
          <Col>
            <SetAddOrderInput addorder={addorder} />
          </Col>
        </Row>
        <Row>
          <Col>
            <SetSolventInput solvent={solvent} />
          </Col>
          <Col>
            <SetConcentrationInput concentration={concentration} />
          </Col>
        </Row>
      </Container>
    </Card>
  );
};

const ReactionNameCard = ({ name }) => {
  const [isLoading, setLoading] = useState(true);
  const [ReactionName, setReactionName] = useState({ name });

  // Must inlcude hook to handle change plus come up with
  // better way of listing names => check with Anthony

  return (
    <Card style={{ width: "8rem" }}>
      <Card.Body>{ReactionName.name}</Card.Body>
    </Card>
  );
};

// Render reactant cards!

// Work up card
const WorkUpCard = ({ workup }) => {
  // Set state
  const [WorkUpType, setWorkUpType] = useState(workup);

  const handleWorkUpChange = (e) => {
    setWorkUpType(e.target.value);
    // Add code to update API as well
  };

  // Render to return
  return (
    <Card style={{ width: "8rem" }}>
      <Card.Img variant="top" src="" alt="Work up image" />
      <Card.Body>
        <Form>
          <Form.Group controlId="exampleForm.ControlSelect1">
            <Form.Control
              as="select"
              onChange={handleWorkUpChange}
              size="sm"
              type="text"
            >
              <option>{WorkUpType.workup}</option>
              <option>DMA</option>
              <option>DCM</option>
            </Form.Control>
          </Form.Group>
        </Form>
      </Card.Body>
    </Card>
  );
};

const TemperatureCard = ({ temperature }) => {
  // Set state
  const [Temperature, setTemperature] = useState({ temperature });

  const handleTemperatureChange = (e) => {
    setTemperature(e.target.value);
    // Add code to update API as well
  };

  // Render to return
  return (
    <Card style={{ width: "8rem" }}>
      <Card.Body>
        <InputGroup
          className="mb-3"
          value={Temperature.temperature}
          onChange={handleTemperatureChange}
          size="sm"
          type="text"
        >
          <FormControl
            aria-describedby="basic-addon1"
            placeholder={Temperature.temperature}
            size="sm"
          />
        </InputGroup>
      </Card.Body>
    </Card>
  );
};

const ReactionCard = ({ reactionid, temperature, workup, name }) => {
  // Use hooks instead of classes
  const [isLoading, setLoading] = useState(true);
  const [Reactants, setReactants] = useState([]);

  useEffect(() => {
    async function fetchData() {
      try {
        const request = await axios.get(`api/reactants?search=${reactionid}`);
        setReactants(request.data);
        setLoading(false);
      } catch (err) {
        console.log(err);
      }
    }
    fetchData();
  }, []);

  if (isLoading) {
    return <div className="App">Loading...</div>;
  }

  return (
    <ListGroup horizontal>
      <ListGroup.Item>
        <ReactionNameCard name={name} />
      </ListGroup.Item>
      {Reactants.map((reactant) => (
        <ListGroup.Item>
          <ReactantCard image={reactant.image} />
        </ListGroup.Item>
      ))}
      <ListGroup.Item>
        <TemperatureCard temperature={temperature} />
      </ListGroup.Item>
      <ListGroup.Item>
        <WorkUpCard workup={workup} />
      </ListGroup.Item>
    </ListGroup>
  );
};

// Up to here with new web app
// const ActionsCard = (() => {
// useSt
// }
const SetDropwise = ({ dropwise, addactionid }) => {
  const [Dropwise, setDropwise] = useState({ dropwise });
  console.log(Dropwise.dropwise === false);

  async function patchDropwise(value) {
    try {
      const response = await axios.patch(`api/addactions/${addactionid}/`, {
        dropwise: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  const handleDropwiseChange = (e) => {
    setDropwise(e.target.value);
    patchDropwise(e.target.value);
  };

  return (
    <InputGroup size="sm" className="mb-3">
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">
          Add dropwise
        </InputGroup.Text>
      </InputGroup.Prepend>
      <Form.Control
        as="select"
        onChange={(event) => handleDropwiseChange(event)}
        size="sm"
        type="text"
      >
        <option>"True"</option>
        {Dropwise.dropwise === true ? (
          <option>"False"</option>
        ) : Dropwise.dropwise === false ? (
          <option>"True"</option>
        ) : (
          <option>"False"</option>
        )}
      </Form.Control>
    </InputGroup>
  );
};

const SetMolEquivalent = ({ molequivalent, addactionid }) => {
  const [TargetMol, setTargetMol] = useState({ molequivalent });

  async function patchMolEquivalent(value) {
    try {
      const response = await axios.patch(`api/addactions/${addactionid}/`, {
        molequivalents: value,
      });
    } catch (error) {
      console.log(error);
    }
  }

  const handleMolEquivalentChange = (e) => {
    const inputMolEquivalent = e.target.value;

    if (!isNaN(inputMolEquivalent)) {
      setTargetMol(e.target.value);
      patchMolEquivalent(Number(e.target.value));
    } else {
      alert("Please input an integer value");
    }
  };

  return (
    <InputGroup size="sm" className="mb-3">
      <InputGroup.Prepend>
        <InputGroup.Text id="inputGroup-sizing-sm">Mol eq.</InputGroup.Text>
      </InputGroup.Prepend>
      <FormControl
        aria-label="Small"
        aria-describedby="inputGroup-sizing-sm"
        placeholder={TargetMol.molequivalent}
        onChange={(event) => handleMolEquivalentChange(event)}
      />
    </InputGroup>
  );
};

const AddActionList = ({ reactionid }) => {
  const [isLoading, setLoading] = useState(true);
  const [AddActions, setAddActions] = useState([]);
  const [Reactants, setReactants] = useState([]);

  useEffect(() => {
    async function fetchData() {
      try {
        const requestadd = await axios.get(
          `api/addactions?search=${reactionid}`
        );
        setAddActions(requestadd.data);
        const requestreact = await axios.get(
          `api/reactants?search=${reactionid}`
        );
        setReactants(requestreact.data);
        setLoading(false);
      } catch (err) {
        console.log(err);
      }
    }
    fetchData();
  }, []);

  if (isLoading) {
    return <div className="App">Loading...</div>;
  }

  async function patchAddNo(addactionid, value) {
    try {
      const response = await axios.patch(`api/addactions/${addactionid}/`, {
        actionno: value,
      });
      console.log(response);
    } catch (error) {
      console.log(error);
    }
  }

  function handleOnDragEnd(result) {
    const additems = Array.from(AddActions);
    const [reorderedadditem] = additems.splice(result.source.index, 1);
    additems.splice(result.destination.index, 0, reorderedadditem);

    const reactitems = Array.from(Reactants);
    const [reorderedreactitems] = reactitems.splice(result.source.index, 1);
    reactitems.splice(result.destination.index, 0, reorderedreactitems);

    setAddActions(additems);
    setReactants(reactitems);

    const addactionid = parseInt(result.draggableId);
    const newactionno = result.destination.index + 1;

    patchAddNo(addactionid, newactionno);
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
            {AddActions.map((addaction, index) => {
              const addactionid = addaction.id.toString();
              return (
                <Draggable
                  key={addaction.id}
                  draggableId={addactionid}
                  index={index}
                >
                  {(provided) => (
                    <ListGroup.Item
                      ref={provided.innerRef}
                      {...provided.draggableProps}
                      {...provided.dragHandleProps}
                    >
                      Step {index + 1} Add
                      <SetMolEquivalent
                        molequivalent={addaction.molequivalents}
                        addactionid={addaction.id}
                      ></SetMolEquivalent>
                      <Image src={Reactants[index].image} fluid />
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

// const ActionList = ({ reactionid }) => {
//   return;
//   <AddActionList reactionid={reactionid}></AddActionList>;
// };

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
    return <div className="App">Loading...</div>;
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
    return <div className="App">Loading...</div>;
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
            <AddActionList reactionid={reaction.id}></AddActionList>
            {/* <ActionList reactionid={reaction.id}></ActionList> */}
          </Accordion.Collapse>
        </Card>
      ))}
    </Accordion>
  );
};

const MethodCard = ({ method }) => {
  return (
    <CardDeck>
      <Card border="light" style={{ width: "18rem" }} key={method.id}>
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
    return <div className="App">Loading...</div>;
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
