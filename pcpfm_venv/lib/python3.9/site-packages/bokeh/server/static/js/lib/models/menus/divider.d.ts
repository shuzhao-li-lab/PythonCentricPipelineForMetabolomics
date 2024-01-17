import { MenuItem, MenuItemView } from "./menu_item";
import type * as p from "../../core/properties";
export declare class DividerView extends MenuItemView {
    model: Divider;
    render(): void;
}
export declare namespace Divider {
    type Attrs = p.AttrsOf<Props>;
    type Props = MenuItem.Props;
}
export interface Divider extends Divider.Attrs {
}
export declare class Divider extends MenuItem {
    properties: Divider.Props;
    __view_type__: DividerView;
    constructor(attrs?: Partial<Divider.Attrs>);
}
//# sourceMappingURL=divider.d.ts.map