import type { PointGeometry } from "../../core/geometry";
import type { FloatArray, ScreenArray } from "../../core/types";
import type { AreaData } from "./area";
import { Area, AreaView } from "./area";
import type { Context2d } from "../../core/util/canvas";
import type { SpatialIndex } from "../../core/util/spatial";
import * as p from "../../core/properties";
import { StepMode } from "../../core/enums";
import { Selection } from "../selections/selection";
export type HAreaStepData = AreaData & {
    _x1: FloatArray;
    _x2: FloatArray;
    _y: FloatArray;
    sx1: ScreenArray;
    sx2: ScreenArray;
    sy: ScreenArray;
};
export interface HAreaStepView extends HAreaStepData {
}
export declare class HAreaStepView extends AreaView {
    model: HAreaStep;
    visuals: HAreaStep.Visuals;
    protected _index_data(index: SpatialIndex): void;
    protected _step_path(ctx: Context2d, mode: StepMode, sx: ScreenArray, sy: ScreenArray, from_i: number, to_i: number): void;
    protected _render(ctx: Context2d, _indices: number[], data?: HAreaStepData): void;
    scenterxy(i: number): [number, number];
    protected _hit_point(geometry: PointGeometry): Selection;
    protected _map_data(): void;
}
export declare namespace HAreaStep {
    type Attrs = p.AttrsOf<Props>;
    type Props = Area.Props & {
        x1: p.CoordinateSpec;
        x2: p.CoordinateSpec;
        y: p.CoordinateSpec;
        step_mode: p.Property<StepMode>;
    };
    type Visuals = Area.Visuals;
}
export interface HAreaStep extends HAreaStep.Attrs {
}
export declare class HAreaStep extends Area {
    properties: HAreaStep.Props;
    __view_type__: HAreaStepView;
    constructor(attrs?: Partial<HAreaStep.Attrs>);
}
//# sourceMappingURL=harea_step.d.ts.map