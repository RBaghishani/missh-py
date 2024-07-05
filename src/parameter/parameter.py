# parameter.py

from typing import List, Tuple
from input import PairFiles
from spaced.spaced_qmer import SpacedQmer

class FileParameter:
    def __init__(self):
        self.input_files = PairFiles()
        self.spaced_seed: List[Tuple[str, SpacedQmer]] = []
        self.num_prev: int = 0

    def init(self, path_file1: str = "", path_file2: str = "") -> bool:
        return self.input_files.init(path_file1, path_file2)

    def add_spaced_qmer(self, name_type: str, spaced: str):
        self.spaced_seed.append((name_type + "/", SpacedQmer(spaced, self.num_prev)))

    def get_input_files(self) -> PairFiles:
        return self.input_files

    def get_v_spaced(self) -> List[Tuple[str, SpacedQmer]]:
        return self.spaced_seed

    def set_num_prev(self, num_prev: int):
        self.num_prev = num_prev
