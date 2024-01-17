import { UIElement, UIElementView } from "../ui/ui_element";
import type { DOMNodeView } from "../dom/dom_node";
import { DOMNode } from "../dom/dom_node";
import type { StyleSheetLike } from "../../core/dom";
import type { IterViews } from "../../core/build_views";
import type * as p from "../../core/properties";
type Button = UIElement;
declare const Button: typeof UIElement;
export declare class DialogView extends UIElementView {
    model: Dialog;
    protected _content: DOMNodeView | UIElementView;
    children(): IterViews;
    stylesheets(): StyleSheetLike[];
    lazy_initialize(): Promise<void>;
    connect_signals(): void;
    remove(): void;
    render(): void;
}
export declare namespace Dialog {
    type Attrs = p.AttrsOf<Props>;
    type Props = UIElement.Props & {
        title: p.Property<string | DOMNode | null>;
        content: p.Property<string | DOMNode | UIElement>;
        buttons: p.Property<Button[]>;
        modal: p.Property<boolean>;
        closable: p.Property<boolean>;
        draggable: p.Property<boolean>;
    };
}
export interface Dialog extends Dialog.Attrs {
}
export declare class Dialog extends UIElement {
    properties: Dialog.Props;
    __view_type__: DialogView;
    constructor(attrs?: Partial<Dialog.Attrs>);
}
export {};
//# sourceMappingURL=dialog.d.ts.map