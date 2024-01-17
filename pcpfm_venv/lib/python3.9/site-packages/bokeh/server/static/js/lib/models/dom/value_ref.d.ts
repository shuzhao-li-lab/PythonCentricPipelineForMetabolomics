import { Placeholder, PlaceholderView } from "./placeholder";
import type { ColumnarDataSource } from "../sources/columnar_data_source";
import type { Index as DataIndex } from "../../core/util/templating";
import type * as p from "../../core/properties";
export declare class ValueRefView extends PlaceholderView {
    model: ValueRef;
    update(source: ColumnarDataSource, i: DataIndex | null, _vars: object): void;
}
export declare namespace ValueRef {
    type Attrs = p.AttrsOf<Props>;
    type Props = Placeholder.Props & {
        field: p.Property<string>;
    };
}
export interface ValueRef extends ValueRef.Attrs {
}
export declare class ValueRef extends Placeholder {
    properties: ValueRef.Props;
    __view_type__: ValueRefView;
    constructor(attrs?: Partial<ValueRef.Attrs>);
}
//# sourceMappingURL=value_ref.d.ts.map