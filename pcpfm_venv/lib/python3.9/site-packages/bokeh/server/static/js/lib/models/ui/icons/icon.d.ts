import { Model } from "../../../model";
import { DOMComponentView } from "../../../core/dom_view";
import type * as p from "../../../core/properties";
export declare abstract class IconView extends DOMComponentView {
    model: Icon;
    connect_signals(): void;
}
export declare namespace Icon {
    type Attrs = p.AttrsOf<Props>;
    type Props = Model.Props & {
        size: p.Property<number | string>;
    };
}
export interface Icon extends Icon.Attrs {
}
export declare abstract class Icon extends Model {
    properties: Icon.Props;
    __view_type__: IconView;
    constructor(attrs?: Partial<Icon.Attrs>);
}
//# sourceMappingURL=icon.d.ts.map