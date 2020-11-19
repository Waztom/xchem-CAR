import React, { Component } from "react";
import ReactDOM from "react-dom";

import Header from "./Layout/Header";
import Body from "./Body/Body";

class App extends Component {
  render() {
    return (
      <React.Fragment>
        <Header />;
        <Body />
      </React.Fragment>
    );
  }
}

ReactDOM.render(<App />, document.getElementById("app"));
