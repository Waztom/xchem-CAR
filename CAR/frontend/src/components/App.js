import React, { useState } from "react";
import ReactDOM from "react-dom";

import Header from "./Layout/Header";
import Body from "./Body/Body";

const App = () => {
  const [ProjectID, setProjectID] = useState(0);

  function changeProject(projectid) {
    setProjectID(projectid);
  }

  function deleteProject() {
    setProjectID(0);
  }

  return (
    <React.Fragment>
      <Header changeProject={changeProject} deleteProject={deleteProject} />
      <Body ProjectID={ProjectID} key={ProjectID} />
    </React.Fragment>
  );
};

ReactDOM.render(<App />, document.getElementById("app"));
