from pathlib import Path
import yaml

from lbm_wetting.utils import pydantic_schemas as schemas


class InputFileParser:
    def __init__(self, file: Path) -> None:
        self.file = file
        with open(self.file) as f:
            self.inputs = yaml.full_load(f)

    def parse(self):
        self.inputs = schemas.Config(**self.inputs).model_dump()
        return self.inputs
