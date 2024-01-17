import { GestureTool, GestureToolView } from "./gesture_tool";
import { DataRenderer } from "../../renderers/data_renderer";
import type * as p from "../../../core/properties";
import type { PinchEvent, ScrollEvent } from "../../../core/ui_events";
import { Dimensions } from "../../../core/enums";
declare const ZoomTogether: import("../../../core/kinds").Kinds.Enum<"none" | "all" | "cross">;
type ZoomTogether = typeof ZoomTogether["__type__"];
declare const Renderers: import("../../../core/kinds").Kinds.Or<[DataRenderer[], "auto"]>;
type Renderers = typeof Renderers["__type__"];
export declare class WheelZoomToolView extends GestureToolView {
    model: WheelZoomTool;
    _pinch(ev: PinchEvent): void;
    _scroll(ev: ScrollEvent): void;
}
export declare namespace WheelZoomTool {
    type Attrs = p.AttrsOf<Props>;
    type Props = GestureTool.Props & {
        dimensions: p.Property<Dimensions>;
        renderers: p.Property<Renderers>;
        level: p.Property<number>;
        maintain_focus: p.Property<boolean>;
        zoom_on_axis: p.Property<boolean>;
        zoom_together: p.Property<ZoomTogether>;
        speed: p.Property<number>;
    };
}
export interface WheelZoomTool extends WheelZoomTool.Attrs {
}
export declare class WheelZoomTool extends GestureTool {
    properties: WheelZoomTool.Props;
    __view_type__: WheelZoomToolView;
    constructor(attrs?: Partial<WheelZoomTool.Attrs>);
    tool_name: string;
    tool_icon: string;
    event_type: "pinch" | "scroll";
    default_order: number;
    get tooltip(): string;
}
export {};
//# sourceMappingURL=wheel_zoom_tool.d.ts.map