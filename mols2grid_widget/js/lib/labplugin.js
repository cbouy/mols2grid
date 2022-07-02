var plugin = require('./index');
var base = require('@jupyter-widgets/base');

module.exports = {
  id: 'mols2grid_widget:plugin',
  requires: [base.IJupyterWidgetRegistry],
  activate: function(app, widgets) {
      widgets.registerWidget({
          name: 'mols2grid_widget',
          version: plugin.version,
          exports: plugin
      });
  },
  autoStart: true
};

