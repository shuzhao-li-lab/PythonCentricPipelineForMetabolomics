import { Control, ControlView } from "./control";
import type { TooltipView } from "../ui/tooltip";
import { Tooltip } from "../ui/tooltip";
import type { IterViews } from "../../core/build_views";
import type { StyleSheetLike } from "../../core/dom";
import type * as p from "../../core/properties";
export type HTMLInputElementLike = HTMLInputElement | HTMLTextAreaElement | HTMLSelectElement;
export declare abstract class InputWidgetView extends ControlView {
    model: InputWidget;
    protected description: TooltipView | null;
    protected input_el: HTMLInputElementLike;
    protected label_el: HTMLLabelElement;
    desc_el: HTMLElement | null;
    protected group_el: HTMLElement;
    controls(): Generator<HTMLInputElementLike, void, unknown>;
    children(): IterViews;
    lazy_initialize(): Promise<void>;
    remove(): void;
    connect_signals(): void;
    stylesheets(): StyleSheetLike[];
    render(): void;
    change_input(): void;
}
export declare namespace InputWidget {
    type Attrs = p.AttrsOf<Props>;
    type Props = Control.Props & {
        title: p.Property<string>;
        description: p.Property<string | Tooltip | null>;
    };
}
export interface InputWidget extends InputWidget.Attrs {
}
export declare abstract class InputWidget extends Control {
    properties: InputWidget.Props;
    __view_type__: InputWidgetView;
    constructor(attrs?: Partial<InputWidget.Attrs>);
}
//# sourceMappingURL=input_widget.d.ts.map