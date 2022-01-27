import React from 'react';
import Alert from 'react-bootstrap/Alert';

function MolAlert({ show }) {
  return (
    <Alert show={show} variant="warning">
      Structure unchanged - invalid input
    </Alert>
  );
}

export default MolAlert;
