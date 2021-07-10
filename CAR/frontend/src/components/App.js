import React, { useState } from "react";
import ReactDOM from "react-dom";

import Header from "./Layout/Header";
import Body from "./Body/Body";
import ProtocolBody from "./Body/ProtocolBody";

const App = () => {
  const [ProjectID, setProjectID] = useState(0);

  function changeProject(projectid) {
    setProjectID(projectid);
  }

  function deleteProject() {
    setProjectID(0);
  }

  function generateProtocol() {
    setProjectID(0);
    // Need to figure out how to setshow Body vs Prot generator -> setShow in state did not work...
  }

  return (
    <React.Fragment>
      <Header
        changeProject={changeProject}
        deleteProject={deleteProject}
        generateProtocol={generateProtocol}
      />
      <Body ProjectID={ProjectID} key={ProjectID} />
      {/* <ProtocolBody ProjectID={ProjectID} key={ProjectID + 1} /> */}
    </React.Fragment>
  );
};

ReactDOM.render(<App />, document.getElementById("app"));
