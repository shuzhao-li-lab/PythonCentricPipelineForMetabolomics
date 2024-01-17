import { Annotation, AnnotationView } from "./annotation";
import { Dimensional } from "./dimensional";
import { Range } from "../ranges/range";
import { Align, Anchor, Orientation, Location } from "../../core/enums";
import type * as visuals from "../../core/visuals";
import * as mixins from "../../core/property_mixins";
import type * as p from "../../core/properties";
import { BBox } from "../../core/util/bbox";
import type { Context2d } from "../../core/util/canvas";
import type { Size, Layoutable } from "../../core/layout";
import { TextLayout } from "../../core/layout";
import { Grid } from "../../core/layout/grid";
import type { ContinuousAxis, ContinuousAxisView } from "../axes/continuous_axis";
import { Ticker } from "../tickers/ticker";
import type { Scale } from "../scales/scale";
declare const LengthSizing: import("../../core/kinds").Kinds.Enum<"adaptive" | "exact">;
type LengthSizing = typeof LengthSizing["__type__"];
export declare class ScaleBarView extends AnnotationView {
    model: ScaleBar;
    visuals: ScaleBar.Visuals;
    bbox: BBox;
    protected label_layout: TextLayout;
    protected title_layout: TextLayout;
    protected axis_layout: Layoutable;
    protected box_layout: Grid;
    protected axis: ContinuousAxis;
    protected axis_view: ContinuousAxisView;
    protected axis_scale: Scale;
    protected cross_scale: Scale;
    protected range: Range;
    protected _get_size(): Size;
    initialize(): void;
    lazy_initialize(): Promise<void>;
    remove(): void;
    connect_signals(): void;
    update_layout(): void;
    update_geometry(): void;
    protected get horizontal(): boolean;
    protected text_layout(args: {
        text: string;
        location: Location;
        align: Align;
        visuals: visuals.Text;
    }): TextLayout;
    compute_geometry(): void;
    protected _draw_box(ctx: Context2d): void;
    protected _draw_axis(_ctx: Context2d): void;
    protected _draw_text(ctx: Context2d, layout: TextLayout, location: Location): void;
    protected _draw_label(ctx: Context2d): void;
    protected _draw_title(ctx: Context2d): void;
    protected _render(): void;
}
export declare namespace ScaleBar {
    type Attrs = p.AttrsOf<Props>;
    type Props = Annotation.Props & {
        range: p.Property<Range | "auto">;
        unit: p.Property<string>;
        dimensional: p.Property<Dimensional>;
        orientation: p.Property<Orientation>;
        bar_length: p.Property<number>;
        length_sizing: p.Property<LengthSizing>;
        location: p.Property<Anchor>;
        label: p.Property<string>;
        label_align: p.Property<Align>;
        label_location: p.Property<Location>;
        label_standoff: p.Property<number>;
        title: p.Property<string>;
        title_align: p.Property<Align>;
        title_location: p.Property<Location>;
        title_standoff: p.Property<number>;
        margin: p.Property<number>;
        padding: p.Property<number>;
        ticker: p.Property<Ticker>;
    } & Mixins;
    type Mixins = mixins.BarLine & mixins.LabelText & mixins.TitleText & mixins.BorderLine & mixins.BackgroundFill & mixins.BackgroundHatch;
    type Visuals = Annotation.Visuals & {
        bar_line: visuals.Line;
        label_text: visuals.Text;
        title_text: visuals.Text;
        border_line: visuals.Line;
        background_fill: visuals.Fill;
        background_hatch: visuals.Hatch;
    };
}
export interface ScaleBar extends ScaleBar.Attrs {
}
export declare class ScaleBar extends Annotation {
    properties: ScaleBar.Props;
    __view_type__: ScaleBarView;
    constructor(attrs?: Partial<ScaleBar.Attrs>);
}
export {};
//# sourceMappingURL=scale_bar.d.ts.map