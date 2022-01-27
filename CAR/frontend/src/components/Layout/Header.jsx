import React, { useState, useEffect } from 'react';
import axios from 'axios';

import Navbar from 'react-bootstrap/Navbar';
import Nav from 'react-bootstrap/Nav';
import Form from 'react-bootstrap/Form';
import FormControl from 'react-bootstrap/FormControl';
import Button from 'react-bootstrap/Button';
import NavDropdown from 'react-bootstrap/NavDropdown';
import Modal from 'react-bootstrap/Modal';

import { Trash } from 'react-bootstrap-icons';
import { HandThumbsUpFill } from 'react-bootstrap-icons';

function DeleteModal({
  show,
  handleDeleteClose,
  handleDelete,
  id,
  projectName,
}) {
  function handleConfirm() {
    handleDelete(id);
    handleDeleteClose();
  }

  return (
    <Modal show={show} onHide={handleDeleteClose} animation={false}>
      <Modal.Header closeButton>
        <Modal.Title>
          <Trash></Trash>Delete Project
        </Modal.Title>
      </Modal.Header>
      <Modal.Body>Warning this will delete project - {projectName}</Modal.Body>
      <Modal.Footer>
        <Button variant="secondary" onClick={() => handleDeleteClose()}>
          Close
        </Button>
        <Button variant="warning" onClick={() => handleConfirm()}>
          Confirm delete
        </Button>
      </Modal.Footer>
    </Modal>
  );
}

function CreateModal({
  show,
  handleCreateClose,
  handleCreate,
  projectName,
  id,
}) {
  function handleConfirm(id) {
    handleCreateClose();
    handleCreate(id);
  }
  return (
    <Modal show={show} onHide={handleCreateClose} animation={false}>
      <Modal.Header closeButton>
        <Modal.Title>
          <HandThumbsUpFill></HandThumbsUpFill>Generate OT Protocol
        </Modal.Title>
      </Modal.Header>
      <Modal.Body>
        Would you like to create auto-generate protocols for project:{' '}
        {projectName}?
      </Modal.Body>
      <Modal.Footer>
        <Button variant="secondary" onClick={() => handleCreateClose()}>
          Close
        </Button>
        <Button variant="warning" onClick={() => handleConfirm(id)}>
          Create Protocol
        </Button>
      </Modal.Footer>
    </Modal>
  );
}

const Header = ({ changeProject, deleteProject, generateProtocol }) => {
  const [Projects, setProjects] = useState([]);
  const [projectID, setProjectID] = useState([]);
  const [projectName, setProjectName] = useState([]);
  const [showDelete, setDeleteShow] = useState(false);
  const [showCreate, setCreateShow] = useState(false);

  useEffect(() => {
    async function fetchData() {
      const request = await axios.get('api/projects/');
      setProjects(request.data);
    }
    fetchData();
  }, []);

  async function deleteData(projectid) {
    try {
      const response = await axios.delete(`api/projects/${projectid}`);
    } catch (error) {
      console.log(error);
    }
  }

  function handleDelete(id) {
    deleteData(id);
    const newList = Projects.filter((item) => item.id !== id);
    setProjects(newList);
    deleteProject();
  }

  function handleCreate(id) {
    generateProtocol(id);
  }

  function handleChange(projectid) {
    changeProject(projectid);
  }

  function handleCreateShow(projectName) {
    setCreateShow(true);
    setProjectName(projectName);
  }

  function handleDeleteShow(projectName, id) {
    setDeleteShow(true);
    setProjectID(id);
    setProjectName(projectName);
  }

  const handleDeleteClose = () => setDeleteShow(false);
  const handleCreateClose = () => setCreateShow(false);

  return (
    <Navbar bg="light" expand="lg">
      <Navbar.Brand href="/">Chemist Assisted Robotics</Navbar.Brand>
      <Navbar.Toggle aria-controls="basic-navbar-nav" />
      <Navbar.Collapse id="basic-navbar-nav">
        <Nav className="mr-auto">
          <Nav.Link href="/upload/">New Project</Nav.Link>
          <NavDropdown title="Load project" id="collasible-nav-dropdown">
            {Projects.map((project) => (
              <NavDropdown.Item
                key={project.id}
                onClick={() => handleChange(project.id)}
              >
                {project.name}
              </NavDropdown.Item>
            ))}
          </NavDropdown>
          <NavDropdown title="Delete project" id="collasible-nav-dropdown">
            {Projects.map((project) => (
              <React.Fragment>
                <NavDropdown.Item
                  key={project.id + 1}
                  onClick={() => handleDeleteShow(project.name, project.id)}
                >
                  {project.name}
                </NavDropdown.Item>
              </React.Fragment>
            ))}
          </NavDropdown>
          <DeleteModal
            show={showDelete}
            handleDeleteClose={handleDeleteClose}
            handleDelete={handleDelete}
            id={projectID}
            projectName={projectName}
          ></DeleteModal>
          <NavDropdown title="Create OT Protocol" id="collasible-nav-dropdown">
            {Projects.map((project) => (
              <React.Fragment>
                <NavDropdown.Item
                  key={project.id + 2}
                  onClick={() => handleCreateShow(project.name)}
                >
                  {project.name}
                </NavDropdown.Item>
              </React.Fragment>
            ))}
          </NavDropdown>
          <CreateModal
            show={showCreate}
            handleCreateClose={handleCreateClose}
            handleCreate={handleCreate}
            projectName={projectName}
          ></CreateModal>
        </Nav>
        <Form inline>
          <FormControl type="text" placeholder="Search" className="mr-sm-2" />
          <Button variant="outline-success">Search</Button>
        </Form>
      </Navbar.Collapse>
    </Navbar>
  );
};

export default Header;
