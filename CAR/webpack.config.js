module.exports = {
  module: {
    rules: [
      { test: /\.css$/i, use: ["style-loader", "css-loader"] },
      {
        test: /\.js$/,
        exclude: /node_modules/,
        use: {
          loader: "babel-loader",
        },
      },
    ],
  },
  watchOptions: { aggregateTimeout: 200, poll: 1000 },
  devtool: "eval-source-map",
};
