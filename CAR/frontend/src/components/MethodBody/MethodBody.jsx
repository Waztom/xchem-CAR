import React, { useState, useEffect, useRef } from 'react';
import axios from 'axios';
import { DragDropContext, Droppable, Draggable } from 'react-beautiful-dnd';

import Card from 'react-bootstrap/Card';
import ListGroup from 'react-bootstrap/ListGroup';
import Button from 'react-bootstrap/Button';
import Accordion from 'react-bootstrap/Accordion';
import CardDeck from 'react-bootstrap/CardDeck';
import Spinner from 'react-bootstrap/Spinner';
import Container from 'react-bootstrap/Container';
import Row from 'react-bootstrap/Row';
import Col from 'react-bootstrap/Col';
import OverlayTrigger from 'react-bootstrap/OverlayTrigger';
import Image from 'react-bootstrap/Image';

import IBMAddAction from '../Actions/IBMAddAction';
import IBMConcentrateAction from '../Actions/IBMConcentrateAction';
import IBMCollectLayerAction from '../Actions/IBMCollectLayerAction';
import IBMDegasAction from '../Actions/IBMDegasAction';
import IBMDrySolidAction from '../Actions/IBMDrySolidAction';
import IBMDrySolutionAction from '../Actions/IBMDrySolutionAction';
import IBMExtractAction from '../Actions/IBMExtractAction';
import IBMFilterAction from '../Actions/IBMFilterAction';
import IBMMakeSolutionAction from '../Actions/IBMMakeSolutionAction';
import IBMPartitionAction from '../Actions/IBMPartitionAction';
import IBMpHAction from '../Actions/IBMpHAction';
import IBMPhaseSeparationAction from '../Actions/IBMPhaseSeparationAction';
import IBMQuenchAction from '../Actions/IBMQuenchAction';
import IBMRefluxAction from '../Actions/IBMRefluxAction';
import IBMSetTemperatureAction from '../Actions/IBMSetTemperatureAction';
import IBMStirAction from '../Actions/IBMStirAction';
import IBMStoreAction from '../Actions/IBMStoreAction.jsx';
import IBMWaitAction from '../Actions/IBMwaitAction';
import IBMWashAction from '../Actions/IBMWashAction';

import ProductImage from '../Images/ProductImage';
import ReactionImage from '../Images/ReactionImage';

String.prototype.capitalize = function () {
  return this.charAt(0).toUpperCase() + this.slice(1);
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

  function getComponent(action, actionno, updateAction, handleDelete) {
    const components = {
      add: (
        <IBMAddAction
          action={action}
          actionno={actionno}
          updateAction={updateAction}
          handleDelete={handleDelete}
        ></IBMAddAction>
      ),
      'collect-layer': (
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
      'dry-solid': (
        <IBMDrySolidAction
          action={action}
          actionno={actionno}
          updateAction={updateAction}
        ></IBMDrySolidAction>
      ),
      'dry-solution': (
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
      'make-solution': (
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
      'phase-separation': (
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
      'set-temperature': (
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
    console.log(changekey);
    var foundIndex = Actions.findIndex((x) => x.id == actionid);
    Actions[foundIndex][changekey] = changevalue;
  }

  async function deleteAction(actiontype, actionid) {
    try {
      const response = await axios.delete(
        `api/IBM${actiontype}actions/${actionid}/`
      );
    } catch (error) {
      console.log(error);
    }
  }

  function handleDelete(actiontype, actionid) {
    deleteAction(actiontype, actionid);
    const newActionList = Actions.filter((item) => item.id !== actionid);
    setActions(newActionList);
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
                '-' +
                action.id.toString() +
                '-' +
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
                      {getComponent(
                        action,
                        actionno,
                        updateAction,
                        handleDelete
                      )}
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

const ReactionAccordian = ({ methodid }) => {
  // Use hooks instead of classes
  const [isLoading, setLoading] = useState(true);
  const [Reactions, setReactions] = useState([]);
  const [isLoadingActions, setLoadActions] = useState(false);
  const ref = useRef(null);

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

  function loadActions() {
    setLoadActions(true);
  }

  return (
    <Accordion>
      {Reactions.map((reaction) => (
        <Card key={reaction.id}>
          <Card.Header>
            <Accordion.Toggle
              as={Button}
              variant="link"
              eventKey={reaction.id}
              onClick={() => loadActions()}
            >
              <Container fluid="md" className="reaction-container">
                <Row className="reaction-container">
                  <Col md="auto">
                    <OverlayTrigger
                      placement="top"
                      delay={{ show: 0, hide: 0 }}
                      overlay={
                        <Image
                          className="reaction-image"
                          src={reaction.reactionimage}
                          fluid
                        />
                      }
                    >
                      <Button className="reaction-image-button">
                        <ProductImage
                          ref={ref}
                          key={reaction.id}
                          reactionid={reaction.id}
                        ></ProductImage>
                      </Button>
                    </OverlayTrigger>
                  </Col>
                  <Col className="reaction-name">{reaction.reactionclass}</Col>
                </Row>
              </Container>
            </Accordion.Toggle>
          </Card.Header>
          {isLoadingActions && (
            <Accordion.Collapse eventKey={reaction.id}>
              <ActionsList
                key={reaction.id}
                reactionid={reaction.id}
              ></ActionsList>
            </Accordion.Collapse>
          )}
        </Card>
      ))}
    </Accordion>
  );
};

const MethodCard = ({ method }) => {
  return (
    <CardDeck>
      <Card className="synthesiscard" key={method.id}>
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
      <ListGroup className="targetmethods" horizontal>
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
