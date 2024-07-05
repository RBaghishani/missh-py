# spaced_qmer_multi.py

import sys
from typing import List
from spaced_qmer import SpacedQmer

class PosOnes:
    def __init__(self, index_one=sys.maxsize, n_one=0, pos_start=0, n_one_before=0):
        self.index_one = index_one
        self.n_one = n_one
        self.pos_start = pos_start
        self.n_one_before = n_one_before

class PreviousShift:
    def __init__(self):
        self.one_to_change = []
        self.one_to_remove = []
        self.one_to_keep = []
        self.one_exit = 0
        self.shift_min = 0
        self.mask = 0

    def get_size(self) -> int:
        return len(self.one_to_change) + len(self.one_to_remove) + self.one_exit + (0 if not self.one_to_remove else 2)

class PreviusShiftExt(PreviousShift):
    def __init__(self):
        super().__init__()
        self.prev_qmer = sys.maxsize
        self.current_sp_ptr = None
        self.prev_sp_ptr = None

    def is_correct_spaced_previous(self) -> bool:
        return self.prev_qmer != sys.maxsize

class MapUnit:
    def __init__(self):
        self.n_one = []
        self.v_v_pos = []

class SpacedQmerMulti:
    def __init__(self):
        self.v_spaced = []
        self.map_unit = MapUnit()
        self.v_shift_min = []
        self.v_shift_min_rotated = []

    def __getitem__(self, i: int) -> SpacedQmer:
        return self.v_spaced[i]

    def __len__(self) -> int:
        return len(self.v_spaced)

    def init(self, v_spaced: List[SpacedQmer]):
        self.v_spaced = v_spaced
        self.get_shift_min_change()
        self.get_shift_min_change_rotated()

    def get_shift_min_change(self):
        table_shift_change = self.get_table_shift_first_on_second()
        max_size = max(len(shift) for shifts in table_shift_change for shift in shifts)

        self.v_shift_min = [[] for _ in range(len(self.v_spaced))]
        for ss in range(len(self.v_spaced)):
            for i in range(max_size):
                limit_evaluation_previous = ss if i == 0 else len(self.v_spaced)
                for s in range(limit_evaluation_previous):
                    if i < len(table_shift_change[ss][s]):
                        current_shift_min = table_shift_change[ss][s][i]
                        current_size = current_shift_min.get_size()
                        current_shift_min.shift_min = i

                        if len(self.v_shift_min[ss]) <= i:
                            if i > 0:
                                saved_prev_shift_min = self.v_shift_min[ss][i-1]
                                saved_size = saved_prev_shift_min.get_size()
                                if saved_size < current_size and saved_prev_shift_min.is_correct_spaced_previous():
                                    current_shift_min = saved_prev_shift_min
                                    current_size = current_shift_min.get_size()
                            if not (i == 0 and current_size > self.v_spaced[ss].get_weight()):
                                self.v_shift_min[ss].append(current_shift_min)
                        else:
                            saved_shift_min = self.v_shift_min[ss][i]
                            saved_size = saved_shift_min.get_size()
                            if saved_size >= current_size or not saved_shift_min.is_correct_spaced_previous():
                                self.v_shift_min[ss][i] = current_shift_min
                if i == 0 and not self.v_shift_min[ss]:
                    self.v_shift_min[ss].append(PreviusShiftExt())

        for ss in range(len(self.v_shift_min)):
            for i in range(len(self.v_shift_min[ss])):
                self.v_shift_min[ss][i].current_sp_ptr = self.v_spaced[ss]
                if self.v_shift_min[ss][i].is_correct_spaced_previous():
                    self.v_shift_min[ss][i].prev_sp_ptr = self.v_spaced[self.v_shift_min[ss][i].prev_qmer]

    def get_table_shift_first_on_second(self):
        table_shift_change = [[[] for _ in range(len(self.v_spaced))] for _ in range(len(self.v_spaced))]
        for ss in range(len(self.v_spaced)):
            for s in range(len(self.v_spaced)):
                table_shift_change[ss][s] = [PreviusShiftExt() for _ in range(len(self.v_spaced[s].spaced_q))]
                for i in range(len(table_shift_change[ss][s])):
                    table_shift_change[ss][s][i].prev_qmer = s

                init = 0
                for i in range(len(self.v_spaced[s].spaced_q)):
                    find = False
                    for j in range(init, len(self.v_spaced[s].pos_one)):
                        if self.v_spaced[s].pos_one[j] >= i:
                            init = j
                            find = True
                            break
                    if not find:
                        init = len(self.v_spaced[s].pos_one)

                    for j in range(init, len(self.v_spaced[s].pos_one)):
                        if (j-init) < len(self.v_spaced[ss].pos_one):
                            if self.v_spaced[ss].pos_one[j-init] != self.v_spaced[s].pos_one[j] - i:
                                table_shift_change[ss][s][i].one_to_remove.append(j-init)
                                table_shift_change[ss][s][i].one_to_change.append(j-init)

                    table_shift_change[ss][s][i].one_exit = init
        return table_shift_change

    def get_shift_min_change_rotated(self):
        shift_firstid_on_secondid = self.v_shift_min
        size_max = max(len(shift) for shift in shift_firstid_on_secondid)
        self.v_shift_min_rotated = [[PreviusShiftExt() for _ in range(len(self.v_spaced))] for _ in range(size_max)]
        for ss in range(len(shift_firstid_on_secondid)):
            for i in range(size_max):
                if i < len(shift_firstid_on_secondid[ss]):
                    self.v_shift_min_rotated[i][ss] = shift_firstid_on_secondid[ss][i]
                else:
                    self.v_shift_min_rotated[i][ss] = shift_firstid_on_secondid[ss][-1]

def print_shift(shift):
    print("one_to_change= ", shift.one_to_change)
    print("one_to_remove= ", shift.one_to_remove)
    print("one_to_keep= ", shift.one_to_keep)
    print("one_exit= ", shift.one_exit)
    print("shift_min= ", shift.shift_min)
    print("mask= ", format(shift.mask, '042b'))
    print()

def printp(pos):
    print(" ".join(map(str, pos)))
