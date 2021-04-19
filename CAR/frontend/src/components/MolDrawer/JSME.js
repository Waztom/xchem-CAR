import React, { useState } from "react";
import { Jsme } from "jsme-react";
import Button from "react-bootstrap/Button";
import Modal from "react-bootstrap/Modal";

function JSME() {
  const [show, setShow] = useState(false);
  const [Smiles, setSmiles] = useState();

  const handleClose = () => setShow(false);
  const handleShow = () => setShow(true);

  function logSmiles(smiles) {
    console.log(smiles);
    setSmiles(smiles);
    console.log(typeof Smiles);
    window.callback(Smiles);
  }

  return (
    <>
      <Button variant="primary" onClick={handleShow}>
        Launch demo modal
      </Button>

      <Modal show={show} onHide={handleClose}>
        <Modal.Header closeButton>
          <Modal.Title>Modal heading</Modal.Title>
        </Modal.Header>
        <Jsme
          height="300px"
          width="400px"
          options="oldlook,star"
          onChange={logSmiles}
        />
        <Modal.Footer>
          <Button variant="secondary" onClick={handleClose}>
            Close
          </Button>
          <Button variant="primary" onClick={handleClose}>
            Save Changes
          </Button>
        </Modal.Footer>
      </Modal>
    </>
  );
}

export default JSME;
