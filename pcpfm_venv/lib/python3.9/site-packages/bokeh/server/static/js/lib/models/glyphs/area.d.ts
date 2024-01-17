import type { GlyphData } from "./glyph";
import { Glyph, GlyphView } from "./glyph";
import type * as visuals from "../../core/visuals";
import type { Rect } from "../../core/types";
import type { Context2d } from "../../core/util/canvas";
import type * as p from "../../core/properties";
import * as mixins from "../../core/property_mixins";
export type AreaData = GlyphData & p.UniformsOf<Area.Mixins>;
export interface AreaView extends AreaData {
}
export declare abstract class AreaView extends GlyphView {
    model: Area;
    visuals: Area.Visuals;
    draw_legend_for_index(ctx: Context2d, bbox: Rect, _index: number): void;
}
export declare namespace Area {
    type Attrs = p.AttrsOf<Props>;
    type Props = Glyph.Props & Mixins;
    type Mixins = mixins.FillScalar & mixins.HatchScalar;
    type Visuals = Glyph.Visuals & {
        fill: visuals.FillScalar;
        hatch: visuals.HatchScalar;
    };
}
export interface Area extends Area.Attrs {
}
export declare class Area extends Glyph {
    properties: Area.Props;
    __view_type__: AreaView;
    constructor(attrs?: Partial<Area.Attrs>);
}
//# sourceMappingURL=area.d.ts.map