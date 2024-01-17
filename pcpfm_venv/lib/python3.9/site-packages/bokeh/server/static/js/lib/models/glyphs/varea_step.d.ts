import type { PointGeometry } from "../../core/geometry";
import type { FloatArray, ScreenArray } from "../../core/types";
import type { AreaData } from "./area";
import { Area, AreaView } from "./area";
import type { Context2d } from "../../core/util/canvas";
import type { SpatialIndex } from "../../core/util/spatial";
import * as p from "../../core/properties";
import { StepMode } from "../../core/enums";
import { Selection } from "../selections/selection";
export type VAreaStepData = AreaData & {
    _x: FloatArray;
    _y1: FloatArray;
    _y2: FloatArray;
    sx: ScreenArray;
    sy1: ScreenArray;
    sy2: ScreenArray;
};
export interface VAreaStepView extends VAreaStepData {
}
export declare class VAreaStepView extends AreaView {
    model: VAreaStep;
    visuals: VAreaStep.Visuals;
    protected _index_data(index: SpatialIndex): void;
    protected _step_path(ctx: Context2d, mode: StepMode, sx: ScreenArray, sy: ScreenArray, from_i: number, to_i: number): void;
    protected _render(ctx: Context2d, _indices: number[], data?: VAreaStepData): void;
    scenterxy(i: number): [number, number];
    protected _hit_point(geometry: PointGeometry): Selection;
    protected _map_data(): void;
}
export declare namespace VAreaStep {
    type Attrs = p.AttrsOf<Props>;
    type Props = Area.Props & {
        x: p.CoordinateSpec;
        y1: p.CoordinateSpec;
        y2: p.CoordinateSpec;
        step_mode: p.Property<StepMode>;
    };
    type Visuals = Area.Visuals;
}
export interface VAreaStep extends VAreaStep.Attrs {
}
export declare class VAreaStep extends Area {
    properties: VAreaStep.Props;
    __view_type__: VAreaStepView;
    constructor(attrs?: Partial<VAreaStep.Attrs>);
}
//# sourceMappingURL=varea_step.d.ts.map