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

const IBMAddAction = ({ action }) => {
  return <Button>{action.actiontype}</Button>;
};
const IBMCollectLayerAction = ({ action }) => {
  return <Button>{action.actiontype}</Button>;
};
const IBMConcentrateAction = ({ action }) => {
  return <Button>{action.actiontype}</Button>;
};
const IBMDegasAction = ({ action }) => {
  return <Button>{action.actiontype}</Button>;
};
const IBMDrySolidAction = ({ action }) => {
  return <Button>{action.actiontype}</Button>;
};
const IBMDrySolutionAction = ({ action }) => {
  return <Button>{action.actiontype}</Button>;
};
const IBMExtractAction = ({ action }) => {
  return <Button>{action.actiontype}</Button>;
};
const IBMFilterAction = ({ action }) => {
  return <Button>{action.actiontype}</Button>;
};
const IBMMakeSolutionAction = ({ action }) => {
  return <Button>{action.actiontype}</Button>;
};
const IBMPartitionAction = ({ action }) => {
  return <Button>{action.actiontype}</Button>;
};
const IBMpHAction = ({ action }) => {
  return <Button>{action.actiontype}</Button>;
};
const IBMPhaseSeparationAction = ({ action }) => {
  return <Button>{action.actiontype}</Button>;
};
const IBMQuenchAction = ({ action }) => {
  return <Button>{action.actiontype}</Button>;
};
const IBMRefluxAction = ({ action }) => {
  return <Button>{action.actiontype}</Button>;
};
const IBMSetTemperatureAction = ({ action }) => {
  return <Button>{action.actiontype}</Button>;
};
const IBMStirAction = ({ action }) => {
  return <Button>{action.actiontype}</Button>;
};
const IBMStoreAction = ({ action }) => {
  return <Button>{action.actiontype}</Button>;
};
const IBMWaitAction = ({ action }) => {
  return <Button>{action.actiontype}</Button>;
};
const IBMWashAction = ({ action }) => {
  return <Button>{action.actiontype}</Button>;
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
        console.log(sorted);
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
      var actiontype = item.actiontype;
      const newactionno = index + 1;
      patchActionNo(actiontype, item.id, newactionno);
    });
  }

  function getComponent(action) {
    const components = {
      add: <IBMAddAction action={action}></IBMAddAction>,
      "collect-layer": (
        <IBMCollectLayerAction action={action}></IBMCollectLayerAction>
      ),
      concentrate: (
        <IBMConcentrateAction action={action}></IBMConcentrateAction>
      ),
      degas: <IBMDegasAction action={action}></IBMDegasAction>,
      "dry-solid": <IBMDrySolidAction action={action}></IBMDrySolidAction>,
      "dry-solution": (
        <IBMDrySolutionAction action={action}></IBMDrySolutionAction>
      ),
      extract: <IBMExtractAction action={action}></IBMExtractAction>,
      filter: <IBMFilterAction action={action}></IBMFilterAction>,
      "make-solution": (
        <IBMMakeSolutionAction action={action}></IBMMakeSolutionAction>
      ),
      partition: <IBMPartitionAction action={action}></IBMPartitionAction>,
      ph: <IBMpHAction action={action}></IBMpHAction>,
      "phase-separation": (
        <IBMPhaseSeparationAction action={action}></IBMPhaseSeparationAction>
      ),
      quench: <IBMQuenchAction action={action}></IBMQuenchAction>,
      reflux: <IBMRefluxAction action={action}></IBMRefluxAction>,
      "set-temperature": (
        <IBMSetTemperatureAction action={action}></IBMSetTemperatureAction>
      ),
      stir: <IBMStirAction action={action}></IBMStirAction>,
      store: <IBMStoreAction action={action}></IBMStoreAction>,
      wait: <IBMWaitAction action={action}></IBMWaitAction>,
      wash: <IBMWashAction action={action}></IBMWashAction>,
    };

    return components[action.actiontype];
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
              const actionid = action.actiontype + "-" + index.toString();
              return (
                <Draggable key={actionid} draggableId={actionid} index={index}>
                  {(provided) => (
                    <ListGroup.Item
                      ref={provided.innerRef}
                      {...provided.draggableProps}
                      {...provided.dragHandleProps}
                    >
                      {getComponent(action)}
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
            <ActionsList reactionid={reaction.id}></ActionsList>
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
