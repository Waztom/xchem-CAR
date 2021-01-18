import React, { useState, useEffect } from "react";
import axios from "axios";
import ReactDOM from "react-dom";

import Header from "./Layout/Header";
import Body from "./Body/Body";

const App = () => {
  const [isLoading, setLoading] = useState(true);
  const [Projects, setProjects] = useState([]);
  const [ProjectID, setProjectID] = useState(0);

  useEffect(() => {
    async function fetchData() {
      const request = await axios.get("api/projects/");
      setProjects(request.data);
      setLoading(false);
    }
    fetchData();
  }, []);

  if (isLoading) {
    return <div className="App">Loading...</div>;
  }

  function handleChange(projectid) {
    setProjectID(projectid);
  }

  return (
    <React.Fragment>
      <Header Projects={Projects} onChange={handleChange} key={Projects.id} />
      <Body ProjectID={ProjectID} key={ProjectID} />
    </React.Fragment>
  );
};

ReactDOM.render(<App />, document.getElementById("app"));
