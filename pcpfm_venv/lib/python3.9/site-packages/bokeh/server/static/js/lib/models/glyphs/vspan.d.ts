import type { GlyphData } from "./glyph";
import { Glyph, GlyphView } from "./glyph";
import { Selection } from "../selections/selection";
import { LineVector } from "../../core/property_mixins";
import type { PointGeometry, SpanGeometry, RectGeometry } from "../../core/geometry";
import type { FloatArray, ScreenArray, Rect } from "../../core/types";
import type * as visuals from "../../core/visuals";
import type { Context2d } from "../../core/util/canvas";
import type { SpatialIndex } from "../../core/util/spatial";
import * as p from "../../core/properties";
export type VSpanData = GlyphData & p.UniformsOf<VSpan.Mixins> & {
    _x: FloatArray;
    sx: ScreenArray;
    max_line_width: number;
};
export interface VSpanView extends VSpanData {
}
export declare class VSpanView extends GlyphView {
    model: VSpan;
    visuals: VSpan.Visuals;
    after_visuals(): void;
    protected _index_data(index: SpatialIndex): void;
    protected _bounds(bounds: Rect): Rect;
    protected _map_data(): void;
    scenterxy(i: number): [number, number];
    protected _render(ctx: Context2d, indices: number[], data?: VSpanData): void;
    protected _get_candidates(sx0: number, sx1?: number): Iterable<number>;
    protected _find_spans(candidates: Iterable<number>, fn: (sx: number, line_width: number) => boolean): number[];
    protected _hit_point(geometry: PointGeometry): Selection;
    protected _hit_span(geometry: SpanGeometry): Selection;
    protected _hit_rect(geometry: RectGeometry): Selection;
    draw_legend_for_index(ctx: Context2d, bbox: Rect, index: number): void;
}
export declare namespace VSpan {
    type Attrs = p.AttrsOf<Props>;
    type Props = Glyph.Props & {
        x: p.CoordinateSpec;
    } & Mixins;
    type Mixins = LineVector;
    type Visuals = Glyph.Visuals & {
        line: visuals.LineVector;
    };
}
export interface VSpan extends VSpan.Attrs {
}
export declare class VSpan extends Glyph {
    properties: VSpan.Props;
    __view_type__: VSpanView;
    constructor(attrs?: Partial<VSpan.Attrs>);
}
//# sourceMappingURL=vspan.d.ts.map