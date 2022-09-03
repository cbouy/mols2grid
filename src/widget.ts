// Copyright (c) Cedric Bouysset
// Distributed under the terms of the Modified BSD License.

import {
  DOMWidgetModel,
  DOMWidgetView,
  ISerializers,
} from '@jupyter-widgets/base';

import { MODULE_NAME, MODULE_VERSION } from './version';

// Import the CSS
import '../css/widget.css';

export class MolGridModel extends DOMWidgetModel {
  defaults() {
    return {
      ...super.defaults(),
      _model_name: MolGridModel.model_name,
      _model_module: MolGridModel.model_module,
      _model_module_version: MolGridModel.model_module_version,
      _view_name: MolGridModel.view_name,
      _view_module: MolGridModel.view_module,
      _view_module_version: MolGridModel.view_module_version,
      grid_id: "default",
      selection: "{}",
      callback_kwargs: "{}",
      filter_mask: [],
    };
  }

  static serializers: ISerializers = {
    ...DOMWidgetModel.serializers,
    // Add any extra serializers here
  };

  static model_name = 'MolGridModel';
  static model_module = MODULE_NAME;
  static model_module_version = MODULE_VERSION;
  static view_name = 'MolGridView'; // Set to null if no view
  static view_module = MODULE_NAME; // Set to null if no view
  static view_module_version = MODULE_VERSION;
}

export class MolGridView extends DOMWidgetView {
  render() {
    this.el.classList.add('mols2grid-widget');
    let grid_id: string = this.model.get('grid_id');
    let name: string = "_MOLS2GRID_" + grid_id;
    (<any>window)[name] = this.model;
    this.model.on('change:filter_mask', this._trigger_filtering, this);
  }

  private _trigger_filtering() {
    let grid_id: string = this.model.get('grid_id');
    let listObj: any = undefined;
    if (typeof (<any>window).mols2grid_lists !== "undefined") {
      listObj = (<any>window).mols2grid_lists[grid_id];
    } else if (typeof (<any>window).parent.mols2grid_lists !== "undefined") {
      listObj = (<any>window).parent.mols2grid_lists[grid_id];
    } else {
      return;
    }
    let filter_mask = this.model.get("filter_mask");
    if (filter_mask !== []) {
      listObj.filter(function (item: any) {
          return filter_mask[item.values()["mols2grid-id"]];
      });
    }
  }
}
