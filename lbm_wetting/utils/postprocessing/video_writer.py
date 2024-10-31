from pathlib import Path

import imageio.v3
import numpy as np
from PIL import Image


class VideoWriter:
    def __init__(self, file: Path | str, fps: int = 30, upsample_rate: int = 0) -> None:
        self.converter = None
        self.image_stack = None
        self.file = None
        self.fps = fps
        self.format = file.suffix
        self.upsample_rate = upsample_rate

        file = Path(file)
        if file.exists():
            raise FileExistsError(f"File {file} already exists")
        file.parent.mkdir(parents=True, exist_ok=True)

        if self.format == ".mp4":
            self.file = imageio.v3.imopen(file, "w", plugin="pyav")
            self.file.init_video_stream("libx264", fps=self.fps)
        elif self.format == ".gif":
            self.image_stack = []
            self.file = file
        else:
            raise ValueError("File format not supported. Use .mp4 or .gif")

    def write(self, frame: np.ndarray, base_image: np.ndarray | None = None) -> None:
        if self.converter is not None:
            if base_image is not None:
                frame = self.converter.convert(base_image=base_image, input_map=frame)
            else:
                frame = self.converter.convert_structure(input_map=frame)
        else:
            if frame.ndim == 2:
                # Convert to rgb image
                frame = np.stack([frame] * 3, axis=-1)

        if self.upsample_rate > 0:
            frame = self._upsample(frame, self.upsample_rate)

        if self.format == ".mp4":
            self.file.write_frame(frame)
        elif self.format == ".gif":
            self.image_stack.append(frame)

    def close(self) -> None:
        if self.format == ".gif":
            self._write_gif(self.image_stack)
        elif self.format == ".mp4":
            # close the video file
            self.file.close()

    def _write_gif(self, image_stack: list) -> None:
        # Convert images to pil images
        image_stack: list[Image.Image] = [
            Image.fromarray(image) for image in image_stack
        ]
        image_stack[0].save(
            self.file,
            save_all=True,
            append_images=image_stack[1:],
            duration=1000 / self.fps,
            loop=0,
        )
