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
export type HSpanData = GlyphData & p.UniformsOf<HSpan.Mixins> & {
    _y: FloatArray;
    sy: ScreenArray;
    max_line_width: number;
};
export interface HSpanView extends HSpanData {
}
export declare class HSpanView extends GlyphView {
    model: HSpan;
    visuals: HSpan.Visuals;
    after_visuals(): void;
    protected _index_data(index: SpatialIndex): void;
    protected _bounds(bounds: Rect): Rect;
    protected _map_data(): void;
    scenterxy(i: number): [number, number];
    protected _render(ctx: Context2d, indices: number[], data?: HSpanData): void;
    protected _get_candidates(sy0: number, sy1?: number): Iterable<number>;
    protected _find_spans(candidates: Iterable<number>, fn: (sy: number, line_width: number) => boolean): number[];
    protected _hit_point(geometry: PointGeometry): Selection;
    protected _hit_span(geometry: SpanGeometry): Selection;
    protected _hit_rect(geometry: RectGeometry): Selection;
    draw_legend_for_index(ctx: Context2d, bbox: Rect, index: number): void;
}
export declare namespace HSpan {
    type Attrs = p.AttrsOf<Props>;
    type Props = Glyph.Props & {
        y: p.CoordinateSpec;
    } & Mixins;
    type Mixins = LineVector;
    type Visuals = Glyph.Visuals & {
        line: visuals.LineVector;
    };
}
export interface HSpan extends HSpan.Attrs {
}
export declare class HSpan extends Glyph {
    properties: HSpan.Props;
    __view_type__: HSpanView;
    constructor(attrs?: Partial<HSpan.Attrs>);
}
//# sourceMappingURL=hspan.d.ts.map