import type { GlyphData } from "./glyph";
import { Glyph, GlyphView } from "./glyph";
import { Selection } from "../selections/selection";
import { LineVector, FillVector, HatchVector } from "../../core/property_mixins";
import type { PointGeometry, SpanGeometry, RectGeometry } from "../../core/geometry";
import type { FloatArray, Rect } from "../../core/types";
import { ScreenArray } from "../../core/types";
import type * as visuals from "../../core/visuals";
import type { Context2d } from "../../core/util/canvas";
import type { SpatialIndex } from "../../core/util/spatial";
import * as p from "../../core/properties";
export type VStripData = GlyphData & p.UniformsOf<VStrip.Mixins> & {
    _x0: FloatArray;
    _x1: FloatArray;
    sx0: ScreenArray;
    sx1: ScreenArray;
    max_width: number;
};
export interface VStripView extends VStripData {
}
export declare class VStripView extends GlyphView {
    model: VStrip;
    visuals: VStrip.Visuals;
    lazy_initialize(): Promise<void>;
    get sleft(): ScreenArray;
    get sright(): ScreenArray;
    get stop(): ScreenArray;
    get sbottom(): ScreenArray;
    protected _set_data(indices: number[] | null): void;
    protected _index_data(index: SpatialIndex): void;
    protected _bounds(bounds: Rect): Rect;
    protected _map_data(): void;
    scenterxy(i: number): [number, number];
    protected _render(ctx: Context2d, indices: number[], data?: VStripData): void;
    protected _get_candidates(sx0: number, sx1?: number): Iterable<number>;
    protected _find_strips(candidates: Iterable<number>, fn: (sx0: number, sx1: number) => boolean): number[];
    protected _hit_point(geometry: PointGeometry): Selection;
    protected _hit_span(geometry: SpanGeometry): Selection;
    protected _hit_rect(geometry: RectGeometry): Selection;
    draw_legend_for_index(ctx: Context2d, bbox: Rect, index: number): void;
}
export declare namespace VStrip {
    type Attrs = p.AttrsOf<Props>;
    type Props = Glyph.Props & {
        x0: p.CoordinateSpec;
        x1: p.CoordinateSpec;
    } & Mixins;
    type Mixins = LineVector & FillVector & HatchVector;
    type Visuals = Glyph.Visuals & {
        line: visuals.LineVector;
        fill: visuals.FillVector;
        hatch: visuals.HatchVector;
    };
}
export interface VStrip extends VStrip.Attrs {
}
export declare class VStrip extends Glyph {
    properties: VStrip.Props;
    __view_type__: VStripView;
    constructor(attrs?: Partial<VStrip.Attrs>);
}
//# sourceMappingURL=vstrip.d.ts.map