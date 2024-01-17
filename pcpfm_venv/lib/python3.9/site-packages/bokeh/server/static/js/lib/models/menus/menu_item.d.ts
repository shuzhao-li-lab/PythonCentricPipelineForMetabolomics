import { UIElement, UIElementView } from "../ui/ui_element";
import type { StyleSheetLike } from "../../core/dom";
import type * as p from "../../core/properties";
export declare abstract class MenuItemView extends UIElementView {
    model: MenuItem;
    stylesheets(): StyleSheetLike[];
}
export declare namespace MenuItem {
    type Attrs = p.AttrsOf<Props>;
    type Props = UIElement.Props;
}
export interface MenuItem extends MenuItem.Attrs {
}
export declare abstract class MenuItem extends UIElement {
    properties: MenuItem.Props;
    __view_type__: MenuItemView;
    constructor(attrs?: Partial<MenuItem.Attrs>);
}
//# sourceMappingURL=menu_item.d.ts.map