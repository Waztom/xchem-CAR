import React from "react";

import Navbar from "react-bootstrap/Navbar";
import Nav from "react-bootstrap/Nav";
import Form from "react-bootstrap/Form";
import FormControl from "react-bootstrap/FormControl";
import Button from "react-bootstrap/Button";
import Dropdown from "react-bootstrap/Dropdown";

const Header = ({ Projects, onChange }) => {
  function handleChange(projectid) {
    onChange(projectid);
  }

  return (
    <Navbar bg="light" expand="lg">
      <Navbar.Brand href="#home">Chemist Assisted Robotics</Navbar.Brand>
      <Navbar.Toggle aria-controls="basic-navbar-nav" />
      <Navbar.Collapse id="basic-navbar-nav">
        <Nav className="mr-auto">
          <Nav.Link href="/upload/">New Project</Nav.Link>
          <Dropdown>
            <Dropdown.Toggle variant="success" id="dropdown-basic">
              Load Project
            </Dropdown.Toggle>
            <Dropdown.Menu>
              {Projects.map((project) => (
                <Dropdown.Item
                  onClick={() => handleChange(project.id)}
                  key={project.id}
                >
                  {project.name}
                </Dropdown.Item>
              ))}
            </Dropdown.Menu>
          </Dropdown>
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
