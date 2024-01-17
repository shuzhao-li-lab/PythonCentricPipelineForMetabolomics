import type { HasProps } from "./has_props";
import type { Property } from "./properties";
import type { Slot, ISignalable } from "./signaling";
import { Signal0, Signal } from "./signaling";
import type { BBox } from "./util/bbox";
export type ViewOf<T extends HasProps> = T["__view_type__"];
export type SerializableState = {
    type: string;
    bbox?: BBox;
    children?: SerializableState[];
};
export declare namespace View {
    type Options = {
        model: HasProps;
        parent: View | null;
        owner?: ViewManager;
    };
}
export type IterViews<T extends View = View> = Generator<T, void, undefined>;
export declare class View implements ISignalable {
    readonly removed: Signal0<this>;
    readonly model: HasProps;
    readonly parent: View | null;
    readonly root: View;
    readonly owner: ViewManager;
    protected _ready: Promise<void>;
    get ready(): Promise<void>;
    children(): IterViews;
    protected _has_finished: boolean;
    mark_finished(): void;
    connect<Args, Sender extends object>(signal: Signal<Args, Sender>, slot: Slot<Args, Sender>): boolean;
    disconnect<Args, Sender extends object>(signal: Signal<Args, Sender>, slot: Slot<Args, Sender>): boolean;
    constructor(options: View.Options);
    initialize(): void;
    lazy_initialize(): Promise<void>;
    protected _removed: boolean;
    remove(): void;
    toString(): string;
    serializable_state(): SerializableState;
    get is_root(): boolean;
    has_finished(): boolean;
    get is_idle(): boolean;
    connect_signals(): void;
    disconnect_signals(): void;
    on_change(properties: Property<unknown> | Property<unknown>[], fn: () => void): void;
    cursor(_sx: number, _sy: number): string | null;
    on_hit?(sx: number, sy: number): boolean;
    private _idle_notified;
    notify_finished(): void;
}
export declare class ViewManager {
    protected drop?: ((view: View) => void) | undefined;
    protected readonly _roots: Set<View>;
    constructor(roots?: Iterable<View>, drop?: ((view: View) => void) | undefined);
    toString(): string;
    get<T extends HasProps>(model: T): ViewOf<T> | null;
    get_by_id(id: string): ViewOf<HasProps> | null;
    add(view: View): void;
    delete(view: View): void;
    clear(): void;
    get roots(): View[];
    [Symbol.iterator](): IterViews;
    views(): IterViews;
    query(fn: (view: View) => boolean): IterViews;
    query_one(fn: (view: View) => boolean): View | null;
    find<T extends HasProps>(model: T): IterViews<ViewOf<T>>;
    find_by_id(id: string): IterViews;
    find_one<T extends HasProps>(model: T): ViewOf<T> | null;
    find_one_by_id(id: string): View | null;
    get_one<T extends HasProps>(model: T): ViewOf<T>;
    get_one_by_id(id: string): View;
    find_all<T extends HasProps>(model: T): ViewOf<T>[];
    find_all_by_id(id: string): View[];
}
//# sourceMappingURL=view.d.ts.map