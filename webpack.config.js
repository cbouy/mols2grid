const path = require('path');
const version = require('./package.json').version;

// Custom webpack rules
const rules = [
  { test: /\.ts$/, loader: 'ts-loader' },
  { test: /\.js$/, loader: 'source-map-loader' },
  { test: /\.css$/, use: ['style-loader', 'css-loader']}
];

// Packages that shouldn't be bundled but loaded at runtime
// 'module' is the magic requirejs dependency used to set the publicPath
const externals = ['@jupyter-widgets/base', 'module'];

const resolve = {
  // Add '.ts' and '.tsx' as resolvable extensions.
  extensions: [".webpack.js", ".web.js", ".ts", ".js"]
};

module.exports = [
  /**
   * Notebook extension
   *
   * This bundle only contains the part of the JavaScript that is run on load of
   * the notebook.
   */
  {
    entry: ['./src/amd-public-path.ts', './src/extension.ts'],
    output: {
      filename: 'index.js',
      path: path.resolve(__dirname, 'mols2grid', 'nbextension'),
      libraryTarget: 'amd',
      publicPath: '', // Set in amd-public-path.js
    },
    module: {
      rules: rules
    },
    devtool: 'source-map',
    externals,
    resolve,
  },

  /**
   * Embeddable mols2grid_widget bundle
   *
   * This bundle is almost identical to the notebook extension bundle. The only
   * difference is in the configuration of the webpack public path for the
   * static assets.
   *
   * The target bundle is always `dist/index.js`, which is the path required by
   * the custom widget embedder.
   */
  {
    entry: ['./src/amd-public-path.ts', './src/index.ts'],
    output: {
        filename: 'index.js',
        path: path.resolve(__dirname, 'dist'),
        libraryTarget: 'amd',
        library: "mols2grid_widget",
        publicPath: '', // Set in amd-public-path.js
    },
    devtool: 'source-map',
    module: {
        rules: rules
    },
    externals,
    resolve,
  },

];
