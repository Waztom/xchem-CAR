import React, { useState, useEffect } from "react";
import axios from "axios";

import Navbar from "react-bootstrap/Navbar";
import Nav from "react-bootstrap/Nav";
import Form from "react-bootstrap/Form";
import FormControl from "react-bootstrap/FormControl";
import Button from "react-bootstrap/Button";
import Dropdown from "react-bootstrap/Dropdown";

import { Trash } from "react-bootstrap-icons";

const Header = ({ onChange }) => {
  const [Projects, setProjects] = useState([]);

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
  }

  function handleChange(projectid) {
    onChange(projectid);
  }

  return (
    <Navbar bg="light" expand="lg">
      <Navbar.Brand href="/">Chemist Assisted Robotics</Navbar.Brand>
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
                <React.Fragment key={project.id}>
                  <Dropdown.Item onClick={() => handleChange(project.id)}>
                    {project.name}
                  </Dropdown.Item>
                  <Button onClick={() => handleDelete(project.id)}>
                    <Trash color="red" size={12}></Trash>
                  </Button>
                </React.Fragment>
              ))}
            </Dropdown.Menu>
          </Dropdown>

          {/* <Nav.Link href="">
            <DeleteProject
              key={ProjectID}
              ProjectID={ProjectID}
              onClick={() => handleChange(0)}
            ></DeleteProject>
          </Nav.Link> */}
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
