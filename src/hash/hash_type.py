from typing import List, Optional

hash_type = int

class Hash_Err:
    def __init__(self):
        self.hash = 0
        self.err_pos: Optional[List[int]] = None

    def reset(self):
        self.hash = 0
        self.err_pos = None

    def sub_pos_err(self, num_pos_to_subtract: int):
        if self.err_pos:
            self.err_pos = [err - num_pos_to_subtract for err in self.err_pos if (err - num_pos_to_subtract) >= 0]

    def add_pos_err(self, num_pos_to_add: int):
        if self.err_pos:
            self.err_pos = [err + num_pos_to_add for err in self.err_pos]

    def sub_pos_err_from(self, num_pos_to_subtract: int, hash_err: 'Hash_Err'):
        if hash_err.err_pos:
            if not self.err_pos:
                self.err_pos = []
            self.err_pos.extend(err - num_pos_to_subtract for err in hash_err.err_pos if (err - num_pos_to_subtract) >= 0)

    def add_pos_err_from(self, num_pos_to_add: int, hash_err: 'Hash_Err'):
        if hash_err.err_pos:
            if not self.err_pos:
                self.err_pos = []
            self.err_pos.extend(err + num_pos_to_add for err in hash_err.err_pos)

    def sort_uniq_err(self):
        if self.err_pos:
            self.err_pos = sorted(set(self.err_pos))

    def is_correct(self) -> bool:
        return not self.err_pos or not self.err_pos

    def create_error(self):
        if not self.err_pos:
            self.err_pos = []

    def push_back_error(self, val: int):
        if not self.err_pos:
            self.err_pos = []
        self.err_pos.append(val)

    def size_error(self) -> int:
        return len(self.err_pos) if self.err_pos else 0

    def __getitem__(self, idx: int) -> int:
        return self.err_pos[idx]

    def __setitem__(self, idx: int, value: int):
        self.err_pos[idx] = value


Hash_Err_V = List[Hash_Err]
Hash_Err_V_V = List[Hash_Err_V]
V_V_Hash_Err = List[Hash_Err_V]
