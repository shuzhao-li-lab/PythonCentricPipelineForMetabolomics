import type { ImageDataBase } from "./image_base";
import { ImageBase, ImageBaseView } from "./image_base";
import { ColorMapper } from "../mappers/color_mapper";
import type { NDArrayType } from "../../core/util/ndarray";
import type * as p from "../../core/properties";
export type ImageData = ImageDataBase;
export interface ImageView extends ImageData {
}
export declare class ImageView extends ImageBaseView {
    model: Image;
    visuals: Image.Visuals;
    load_glglyph(): Promise<typeof import("./webgl/image").ImageGL>;
    connect_signals(): void;
    protected _update_image(): void;
    protected _flat_img_to_buf8(img: NDArrayType<number>): Uint8ClampedArray;
}
export declare namespace Image {
    type Attrs = p.AttrsOf<Props>;
    type Props = ImageBase.Props & {
        color_mapper: p.Property<ColorMapper>;
    };
    type Visuals = ImageBase.Visuals;
}
export interface Image extends Image.Attrs {
}
export declare class Image extends ImageBase {
    properties: Image.Props;
    __view_type__: ImageView;
    constructor(attrs?: Partial<Image.Attrs>);
}
//# sourceMappingURL=image.d.ts.map