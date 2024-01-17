import type { ReglWrapper } from "./regl_wrap";
import type { Float32Buffer } from "./buffer";
import { SXSYGlyphGL } from "./sxsy";
import type { GLMarkerType } from "./types";
import type { CircleView } from "../circle";
export declare class CircleGL extends SXSYGlyphGL {
    readonly glyph: CircleView;
    constructor(regl_wrapper: ReglWrapper, glyph: CircleView);
    get marker_type(): GLMarkerType;
    get size(): Float32Buffer;
    protected _set_data(): void;
    protected _set_once(): void;
}
//# sourceMappingURL=circle.d.ts.map