import type { StyleSheetLike } from "../../core/dom";
import { DropPane } from "../../core/util/panes";
import type * as p from "../../core/properties";
import { InputWidget, InputWidgetView } from "./input_widget";
declare const Item: import("../../core/kinds").Kinds.Tuple<[string, import("../../core/types").Arrayable<import("../../core/types").Color>]>;
type Item = typeof Item["__type__"];
export declare class ColorMapView extends InputWidgetView {
    model: ColorMap;
    input_el: HTMLSelectElement;
    protected _pane: DropPane;
    stylesheets(): StyleSheetLike[];
    connect_signals(): void;
    protected _render_item(item: Item): HTMLDivElement;
    render(): void;
    select(item: Item): void;
    toggle(): void;
    hide(): void;
}
export declare namespace ColorMap {
    type Attrs = p.AttrsOf<Props>;
    type Props = InputWidget.Props & {
        value: p.Property<string>;
        items: p.Property<Item[]>;
        swatch_width: p.Property<number>;
        swatch_height: p.Property<number>;
        ncols: p.Property<number>;
    };
}
export interface ColorMap extends ColorMap.Attrs {
}
export declare class ColorMap extends InputWidget {
    properties: ColorMap.Props;
    __view_type__: ColorMapView;
    constructor(attrs?: Partial<ColorMap.Attrs>);
}
export {};
//# sourceMappingURL=color_map.d.ts.map