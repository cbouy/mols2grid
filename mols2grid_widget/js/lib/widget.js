var widgets = require('@jupyter-widgets/base');
var _ = require('lodash');

var CommWidgetModel = widgets.DOMWidgetModel.extend({
    defaults: _.extend(widgets.DOMWidgetModel.prototype.defaults(), {
        _model_name : 'CommWidgetModel',
        _view_name : 'CommWidgetView',
        _model_module : 'mols2grid_widget',
        _view_module : 'mols2grid_widget',
        _model_module_version : '0.1.0',
        _view_module_version : '0.1.0',
        grid_id : "default",
        selection : "{}",
        callback_args : "{}"
    })
});


var CommWidgetView = widgets.DOMWidgetView.extend({
    // Defines how the widget gets rendered into the DOM
    render: function() {
        let grid_id = this.model.get('grid_id');
        let name = "_MOLS2GRID_" + grid_id;
        window[name] = this.model;
        // Python -> JS
        this.model.on('change:grid_id', this.value_changed, this);
    },

    value_changed: function() {
        
    }
});


module.exports = {
    CommWidgetModel: CommWidgetModel,
    CommWidgetView: CommWidgetView
};
