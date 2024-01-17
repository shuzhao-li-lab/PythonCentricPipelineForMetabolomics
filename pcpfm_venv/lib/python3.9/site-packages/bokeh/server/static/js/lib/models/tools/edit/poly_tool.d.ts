import type * as p from "../../../core/properties";
import type { UIEvent } from "../../../core/ui_events";
import type { MultiLine } from "../../glyphs/multi_line";
import type { Patches } from "../../glyphs/patches";
import type { GlyphRenderer } from "../../renderers/glyph_renderer";
import type { HasXYGlyph } from "./edit_tool";
import { EditTool, EditToolView } from "./edit_tool";
export interface HasPolyGlyph {
    glyph: MultiLine | Patches;
}
export declare abstract class PolyToolView extends EditToolView {
    model: PolyTool;
    _set_vertices(xs: number[] | number, ys: number[] | number): void;
    _hide_vertices(): void;
    _snap_to_vertex(ev: UIEvent, x: number, y: number): [number, number];
}
export declare namespace PolyTool {
    type Attrs = p.AttrsOf<Props>;
    type Props = EditTool.Props & {
        vertex_renderer: p.Property<(GlyphRenderer & HasXYGlyph) | null>;
    };
}
export interface PolyTool extends PolyTool.Attrs {
}
export declare abstract class PolyTool extends EditTool {
    properties: PolyTool.Props;
    __view_type__: PolyToolView;
    constructor(attrs?: Partial<PolyTool.Attrs>);
}
//# sourceMappingURL=poly_tool.d.ts.map