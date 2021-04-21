import React, { useState } from "react";
import axios from "axios";

import Image from "react-bootstrap/Image";
import Container from "react-bootstrap/Container";
import Row from "react-bootstrap/Row";
import Col from "react-bootstrap/Col";
import Button from "react-bootstrap/Button";
import { Trash } from "react-bootstrap-icons";

import SetQuantity from "../SetActionInputs/SetQuantity";
import SetDropwise from "../SetActionInputs/SetDropwise";
import SetAtmosphere from "../SetActionInputs/SetAtmosphere";

import JSMEModal from "../MolDrawer/JSMEModal";
import MolAlert from "../MolDrawer/MolAlert";

import "../styles.css";

const IBMAddAction = ({ action, actionno, updateAction, handleDelete }) => {
  const actiontype = action.actiontype.capitalize();
  const id = action.id;
  const [Show, setShow] = useState(false);
  const [ShowAlert, setShowAlert] = useState(false);
  const [Smiles, setSmiles] = useState(action.materialsmiles);
  const [SVG, setSVG] = useState(action.materialimage);

  function handleShow() {
    setShow(true);
  }

  function handleClose(smiles) {
    setShow(false);
    var mol = window.checkmol(smiles);
    if (mol) {
      setShowAlert(false);
      var origmol = window.checkmol(Smiles);
      var origsmiles = origmol.get_smiles();
      var canonsmiles = mol.get_smiles();
    } else {
      setShowAlert(true);
      console.log("mol invalid");
    }
    if (canonsmiles !== origsmiles) {
      setSmiles(canonsmiles);
      var svg = window.getsvg(mol);
      var blob = new Blob([svg], { type: "image/svg+xml" });
      var url = URL.createObjectURL(blob);
      setSVG(url);
      patchSVG(canonsmiles);
    }
  }

  async function patchSVG(value) {
    try {
      const response = await axios.patch(
        `api/IBM${action.actiontype}actions/${id}/`,
        {
          materialsmiles: value,
          changeimage: 0,
        }
      );
    } catch (error) {
      console.log(error);
    }
    console.log(value);
  }

  return (
    <Container key={actionno}>
      <h5>
        {actionno}. {actiontype}
      </h5>
      <Row>
        <Col>
          <MolAlert show={ShowAlert}></MolAlert>
          <JSMEModal
            show={Show}
            handleClose={handleClose}
            smiles={Smiles}
          ></JSMEModal>
          <Button className="editcompound" onClick={() => handleShow()}>
            <Image width={300} height={100} src={SVG} alt={action.material} />
          </Button>
          <SetQuantity
            action={action}
            updateAction={updateAction}
            name={"material"}
          ></SetQuantity>
          <SetDropwise
            action={action}
            updateAction={updateAction}
          ></SetDropwise>
          <SetAtmosphere
            action={action}
            updateAction={updateAction}
          ></SetAtmosphere>
        </Col>
      </Row>
      <Row>
        <Col>
          <Button
            key={action.id}
            onClick={() => handleDelete(action.actiontype, action.id)}
          >
            <Trash></Trash>
          </Button>
        </Col>
      </Row>
    </Container>
  );
};

export default IBMAddAction;
