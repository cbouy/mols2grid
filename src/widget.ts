import type { AnyModel, RenderProps } from "@anywidget/types";
import "./widget.css";

/* Specifies attributes defined with traitlets in widget.py */
interface WidgetModel {
	grid_id: string,
	selection: string,
	callback_kwargs: string,
	filter_mask: boolean[],
}

function render({ model, el }: RenderProps<WidgetModel>) {
	let grid_id: string = model.get("grid_id");
    let name: string = "_MOLS2GRID_" + grid_id;
    (<any>window)[name] = model;
    model.on("change:filter_mask", () => {
        trigger_filtering(model)
    });
	el.classList.add("mols2grid-anywidget");
}

function trigger_filtering(model: AnyModel<WidgetModel>) {
	let grid_id: string = model.get("grid_id");
    let listObj: any = undefined;
    if (typeof (<any>window).mols2grid_lists !== "undefined") {
      	listObj = (<any>window).mols2grid_lists[grid_id];
    } else if (typeof (<any>window).parent.mols2grid_lists !== "undefined") {
      	listObj = (<any>window).parent.mols2grid_lists[grid_id];
    } else {
      	return;
    }
    let filter_mask = model.get("filter_mask");
    if (filter_mask) {
		listObj.filter(function (item: any) {
			return filter_mask[item.values()["mols2grid-id"]];
		});
    }
}

export default { render };
