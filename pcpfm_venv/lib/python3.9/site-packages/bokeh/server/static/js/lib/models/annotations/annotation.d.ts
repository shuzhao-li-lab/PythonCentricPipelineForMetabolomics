import { Renderer, RendererView } from "../renderers/renderer";
import type { Panel } from "../../core/layout/side_panel";
import type { Size, Layoutable } from "../../core/layout";
import type { SerializableState } from "../../core/view";
import type { BBox } from "../../core/util/bbox";
import type * as p from "../../core/properties";
export declare abstract class AnnotationView extends RendererView {
    model: Annotation;
    layout?: Layoutable;
    panel?: Panel;
    bbox?: BBox;
    update_layout?(): void;
    after_layout?(): void;
    get_size(): Size;
    protected _get_size(): Size;
    connect_signals(): void;
    get needs_clip(): boolean;
    serializable_state(): SerializableState;
}
export declare namespace Annotation {
    type Attrs = p.AttrsOf<Props>;
    type Props = Renderer.Props;
    type Visuals = Renderer.Visuals;
}
export interface Annotation extends Annotation.Attrs {
}
export declare abstract class Annotation extends Renderer {
    properties: Annotation.Props;
    __view_type__: AnnotationView;
    constructor(attrs?: Partial<Annotation.Attrs>);
}
//# sourceMappingURL=annotation.d.ts.map