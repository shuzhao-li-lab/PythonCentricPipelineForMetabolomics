import { Model } from "../../model";
import type * as p from "../../core/properties";
export type BasisItem = {
    factor: number;
    short_name: string;
    long_name: string;
    tex_repr: string;
};
export type Basis = BasisItem[];
export type PreferredValue = {
    new_value: number;
    new_unit: string;
    new_long_unit: string;
    scale_factor: number;
    exact: boolean;
};
export declare namespace Dimensional {
    type Attrs = p.AttrsOf<Props>;
    type Props = Model.Props & {
        ticks: p.Property<number[]>;
        include: p.Property<string[] | null>;
        exclude: p.Property<string[]>;
    };
}
export interface Dimensional extends Dimensional.Attrs {
}
export declare abstract class Dimensional extends Model {
    properties: Dimensional.Props;
    constructor(attrs?: Partial<Dimensional.Attrs>);
    abstract get basis(): Basis;
    compute(value: number, unit: string, exact?: boolean): PreferredValue;
}
export declare namespace Metric {
    type Attrs = p.AttrsOf<Props>;
    type Props = Dimensional.Props;
}
export interface Metric extends Metric.Attrs {
}
export declare abstract class Metric extends Dimensional {
    properties: Metric.Props;
    constructor(attrs?: Partial<Metric.Attrs>);
    build_basis(): Basis;
    protected _basis: Basis | null;
    get basis(): Basis;
    protected abstract _short_name: string;
    protected abstract _long_name: string;
    protected abstract _tex_repr: string;
    protected _basis_template: [short_name: string, factor: number, long_name: string, tex: string | null][];
}
export declare namespace MetricLength {
    type Attrs = p.AttrsOf<Props>;
    type Props = Metric.Props;
}
export interface MetricLength extends MetricLength.Attrs {
}
export declare class MetricLength extends Metric {
    properties: MetricLength.Props;
    constructor(attrs?: Partial<MetricLength.Attrs>);
    protected _short_name: string;
    protected _long_name: string;
    protected _tex_repr: string;
}
//# sourceMappingURL=dimensional.d.ts.map