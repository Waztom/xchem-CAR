import React, { useState, useEffect } from "react";
import axios from "axios";

import Navbar from "react-bootstrap/Navbar";
import Nav from "react-bootstrap/Nav";
import Form from "react-bootstrap/Form";
import FormControl from "react-bootstrap/FormControl";
import Button from "react-bootstrap/Button";
import NavDropdown from "react-bootstrap/NavDropdown";
import Modal from "react-bootstrap/Modal";

import { Trash } from "react-bootstrap-icons";

function DeleteModal({ show, handleClose, handleDelete, id }) {
  function handleConfirm() {
    handleDelete(id);
    handleClose();
  }

  return (
    <Modal show={show} onHide={handleClose} animation={false}>
      <Modal.Header closeButton>
        <Modal.Title>
          <Trash></Trash>Delete Project
        </Modal.Title>
      </Modal.Header>
      <Modal.Body>Warning this will delete your project!</Modal.Body>
      <Modal.Footer>
        <Button variant="secondary" onClick={() => handleClose()}>
          Close
        </Button>
        <Button variant="warning" onClick={() => handleConfirm()}>
          Confirm delete
        </Button>
      </Modal.Footer>
    </Modal>
  );
}

const Header = ({ changeProject, deleteProject }) => {
  const [Projects, setProjects] = useState([]);
  const [show, setShow] = useState(false);

  useEffect(() => {
    async function fetchData() {
      const request = await axios.get("api/projects/");
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

  function handleChange(projectid) {
    changeProject(projectid);
  }

  const handleClose = () => setShow(false);
  const handleShow = () => setShow(true);

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
                <NavDropdown.Item key={project.id} onClick={() => handleShow()}>
                  {project.name}
                </NavDropdown.Item>
                <DeleteModal
                  show={show}
                  handleClose={handleClose}
                  handleDelete={handleDelete}
                  id={project.id}
                ></DeleteModal>
              </React.Fragment>
            ))}
          </NavDropdown>
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
