from pathlib import Path

import imageio.v3
import numpy as np
from PIL import Image


class VideoWriter:
    def __init__(self, file: Path | str, fps: int = 30, upsample_rate: int = 0) -> None:
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

    def write(self, frame: np.ndarray) -> None:
        # convert color of the frame
        colored_frame = np.copy(frame)
        colored_frame[frame == 255] = 0
        colored_frame[frame == 200] = 100
        colored_frame[frame == 0] = 255

        colored_frame = np.stack([colored_frame] * 3, axis=-1, dtype=np.uint8)
        colored_frame[frame == 128] = np.array([0, 86, 158], dtype=np.uint8)

        if self.upsample_rate > 0:
            colored_frame = self._upsample(colored_frame, self.upsample_rate)

        if self.format == ".mp4":
            self.file.write_frame(colored_frame)
        elif self.format == ".gif":
            self.image_stack.append(colored_frame)

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

    def _upsample(self, frame: np.ndarray, upsample_rate: int) -> np.ndarray:
        return np.kron(
            frame, np.ones((upsample_rate, upsample_rate, 1), dtype=np.uint8)
        )
