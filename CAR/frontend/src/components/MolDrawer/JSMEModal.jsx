import React, { useState } from 'react';
import { Jsme } from 'jsme-react';
import Button from 'react-bootstrap/Button';
import Modal from 'react-modal';

import '../styles.css';

// Make sure to bind modal to your appElement (https://reactcommunity.org/react-modal/accessibility/)
// Modal.setAppElement(document.getElementById("root"));

function JSMEModal({ show, smiles, handleClose }) {
  const [Smiles, setSmiles] = useState(smiles);

  function logSmiles(smiles) {
    setSmiles(smiles);
  }

  return (
    <Modal isOpen={show} className="drawmodal" ariaHideApp={false}>
      <Jsme width="600px" height="700px" smiles={smiles} onChange={logSmiles} />
      <Button onClick={() => handleClose(Smiles)}>close</Button>
    </Modal>
  );
}

export default JSMEModal;
