import type * as p from "../../../core/properties";
import type { Line } from "../../glyphs/line";
import type { GlyphRenderer } from "../../renderers/glyph_renderer";
import { EditTool, EditToolView } from "./edit_tool";
export type HasLineGlyph = {
    glyph: Line;
};
export declare abstract class LineToolView extends EditToolView {
    model: LineTool;
    _set_intersection(x: number[] | number, y: number[] | number): void;
    _hide_intersections(): void;
}
export declare namespace LineTool {
    type Attrs = p.AttrsOf<Props>;
    type Props = EditTool.Props & {
        intersection_renderer: p.Property<(GlyphRenderer & HasLineGlyph)>;
    };
}
export interface LineTool extends LineTool.Attrs {
}
export declare abstract class LineTool extends EditTool {
    properties: LineTool.Props;
    __view_type__: LineToolView;
    constructor(attrs?: Partial<LineTool.Attrs>);
}
//# sourceMappingURL=line_tool.d.ts.map