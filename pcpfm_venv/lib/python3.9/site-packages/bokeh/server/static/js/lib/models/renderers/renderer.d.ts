import { View } from "../../core/view";
import * as visuals from "../../core/visuals";
import { RenderLevel } from "../../core/enums";
import type * as p from "../../core/properties";
import { Model } from "../../model";
import type { CanvasLayer } from "../../core/util/canvas";
import type { Plot, PlotView } from "../plots/plot";
import type { CanvasView } from "../canvas/canvas";
import { CoordinateTransform, CoordinateMapping } from "../coordinates/coordinate_mapping";
import type { Node } from "../coordinates/node";
export declare namespace RendererGroup {
    type Attrs = p.AttrsOf<Props>;
    type Props = Model.Props & {
        visible: p.Property<boolean>;
    };
}
export interface RendererGroup extends RendererGroup.Attrs {
}
export declare class RendererGroup extends Model {
    properties: RendererGroup.Props;
    constructor(attrs?: Partial<RendererGroup.Attrs>);
}
export declare abstract class RendererView extends View implements visuals.Renderable {
    model: Renderer;
    visuals: Renderer.Visuals;
    readonly parent: PlotView;
    protected _coordinates?: CoordinateTransform;
    get coordinates(): CoordinateTransform;
    private _custom_coordinates;
    set coordinates(custom_coordinates: CoordinateTransform | null);
    initialize(): void;
    connect_signals(): void;
    protected _initialize_coordinates(): CoordinateTransform;
    get plot_view(): PlotView;
    get plot_model(): Plot;
    get layer(): CanvasLayer;
    get canvas(): CanvasView;
    request_render(): void;
    request_paint(): void;
    request_layout(): void;
    notify_finished(): void;
    notify_finished_after_paint(): void;
    interactive_hit?(sx: number, sy: number): boolean;
    get needs_clip(): boolean;
    get has_webgl(): boolean;
    get displayed(): boolean;
    render(): void;
    protected abstract _render(): void;
    renderer_view<T extends Renderer>(_renderer: T): T["__view_type__"] | undefined;
    /**
     * Geometry setup that doesn't change between paints.
     */
    update_geometry(): void;
    /**
     * Geometry setup that changes between paints.
     */
    compute_geometry(): void;
    /**
     * Compute screen coordinates for a symbolic node.
     */
    resolve_node(node: Node): {
        x: number;
        y: number;
    };
}
export declare namespace Renderer {
    type Attrs = p.AttrsOf<Props>;
    type Props = Model.Props & {
        group: p.Property<RendererGroup | null>;
        level: p.Property<RenderLevel>;
        visible: p.Property<boolean>;
        x_range_name: p.Property<string>;
        y_range_name: p.Property<string>;
        coordinates: p.Property<CoordinateMapping | null>;
        propagate_hover: p.Property<boolean>;
    };
    type Visuals = visuals.Visuals;
}
export interface Renderer extends Renderer.Attrs {
}
export declare abstract class Renderer extends Model {
    properties: Renderer.Props;
    __view_type__: RendererView;
    constructor(attrs?: Partial<Renderer.Attrs>);
}
//# sourceMappingURL=renderer.d.ts.map