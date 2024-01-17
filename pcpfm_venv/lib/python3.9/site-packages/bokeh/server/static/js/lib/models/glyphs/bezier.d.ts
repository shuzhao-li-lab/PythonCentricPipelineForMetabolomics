import { LineVector } from "../../core/property_mixins";
import type * as visuals from "../../core/visuals";
import type { Rect, FloatArray, ScreenArray } from "../../core/types";
import type { SpatialIndex } from "../../core/util/spatial";
import type { Context2d } from "../../core/util/canvas";
import type { GlyphData } from "./glyph";
import { Glyph, GlyphView } from "./glyph";
import * as p from "../../core/properties";
export type BezierData = GlyphData & p.UniformsOf<Bezier.Mixins> & {
    _x0: FloatArray;
    _y0: FloatArray;
    _x1: FloatArray;
    _y1: FloatArray;
    _cx0: FloatArray;
    _cy0: FloatArray;
    _cx1: FloatArray;
    _cy1: FloatArray;
    sx0: ScreenArray;
    sy0: ScreenArray;
    sx1: ScreenArray;
    sy1: ScreenArray;
    scx0: ScreenArray;
    scy0: ScreenArray;
    scx1: ScreenArray;
    scy1: ScreenArray;
};
export interface BezierView extends BezierData {
}
export declare class BezierView extends GlyphView {
    model: Bezier;
    visuals: Bezier.Visuals;
    protected _project_data(): void;
    protected _index_data(index: SpatialIndex): void;
    protected _render(ctx: Context2d, indices: number[], data?: BezierData): void;
    draw_legend_for_index(ctx: Context2d, bbox: Rect, index: number): void;
    scenterxy(): [number, number];
}
export declare namespace Bezier {
    type Attrs = p.AttrsOf<Props>;
    type Props = Glyph.Props & {
        x0: p.CoordinateSpec;
        y0: p.CoordinateSpec;
        x1: p.CoordinateSpec;
        y1: p.CoordinateSpec;
        cx0: p.CoordinateSpec;
        cy0: p.CoordinateSpec;
        cx1: p.CoordinateSpec;
        cy1: p.CoordinateSpec;
    } & Mixins;
    type Mixins = LineVector;
    type Visuals = Glyph.Visuals & {
        line: visuals.LineVector;
    };
}
export interface Bezier extends Bezier.Attrs {
}
export declare class Bezier extends Glyph {
    properties: Bezier.Props;
    __view_type__: BezierView;
    constructor(attrs?: Partial<Bezier.Attrs>);
}
//# sourceMappingURL=bezier.d.ts.map