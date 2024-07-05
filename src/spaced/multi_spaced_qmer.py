from typing import List, NamedTuple

class Position(NamedTuple):
    positions: List[int]

class PreviousShiftMulti(NamedTuple):
    one_to_remove: Position
    one_to_keep: Position
    one_exit: int = 0
    shift_min: int = 0
    seed_num: int = 0
    mask: int = 0

class groupPrevious(NamedTuple):
    prev: List[PreviousShiftMulti]
    not_covered: Position

class SeedInfo(NamedTuple):
    group_previous: List[groupPrevious]
    pos_ones: Position

class MultiSeedInfo(NamedTuple):
    info: List[SeedInfo]

class MultiSeedInfoRow(NamedTuple):
    info: List[SeedInfo]
    transient1: int
    transient2: int

class MultiSpacedQmer:
    def __init__(self, spaced_qmers: List['SpacedQmer']):
        self.spaced_qmers = spaced_qmers
        self.spaced_qs = [None] * self.get_length()
        self.v_pos_ones = [None] * self.get_length()
        self.v_pos_pos_ones = [None] * self.get_length()
        length = spaced_qmers[0].get_q()
        weight = spaced_qmers[0].get_weight()

        for i in range(self.get_length()):
            if length != spaced_qmers[i].get_q() or weight != spaced_qmers[i].get_weight():
                raise ValueError("The seed in the group does not have the same length or weight")

            self.spaced_qs[i] = spaced_qmers[i].to_string()
            k = 0
            for j in range(len(self.spaced_qs[i])):
                if spaced_qmers[i].is_one(j):
                    if self.v_pos_ones[i] is None:
                        self.v_pos_ones[i] = []
                    if self.v_pos_pos_ones[i] is None:
                        self.v_pos_pos_ones[i] = []
                    self.v_pos_ones[i].append(j)
                    self.v_pos_pos_ones[i].append(k)
                    k += 1

        self.set_multi_seed_info_row()
        self.v_pos_pos_ones = [None] * self.get_length()
        for i in range(self.get_length()):
            k = 0
            for j in range(len(self.spaced_qs[i])):
                if spaced_qmers[i].is_one(j):
                    if self.v_pos_pos_ones[i] is None:
                        self.v_pos_pos_ones[i] = []
                    self.v_pos_pos_ones[i].append(k)
                    k += 1
        self.set_multi_seed_info_col()

    def get_length(self):
        return len(self.spaced_qmers)

    def get_multi_seed_info_col(self):
        return self.multi_seed_info_col

    def get_multi_seed_info_row(self):
        return self.multi_seed_info_row

    def reset(self, spaced_qmers: List['SpacedQmer']):
        self.__init__(spaced_qmers)

    def set_multi_seed_info_col(self):
        self.multi_seed_info_col = [None] * self.get_length()

        for j in range(self.get_length()):
            self.multi_seed_info_col[j] = SeedInfo(group_previous=[groupPrevious(prev=[], not_covered=Position([]))], pos_ones=Position(self.v_pos_ones[j]))

        self.process_multi_seed_col(0)

        furthest_pos = 0
        for j in range(self.get_length()):
            for prev in self.multi_seed_info_col[j].group_previous[0].prev:
                furthest_pos = max(furthest_pos, prev.shift_min)
        furthest_pos += 1

        for j in range(self.get_length()):
            self.multi_seed_info_col[j].group_previous = [groupPrevious(prev=[], not_covered=Position([])) for _ in range(furthest_pos)]

        for k in range(1, furthest_pos):
            self.process_multi_seed_col(k)

    def process_multi_seed_col(self, index):
        for y in range(self.get_length()):
            current_shift_group = []
            pos_not_covered_yet = self.v_pos_pos_ones[y]

            done = False
            while not done:
                curr_best = PreviousShiftMulti(Position([]), Position([]), float('inf'))

                max_num_shifts = index - 1 if index != 0 else len(self.spaced_qs[index]) - 1
                for n in range(self.get_length()):
                    starting_point = 0 if n < y else 1
                    for i in range(starting_point, max_num_shifts + 1):
                        starting_init = 1 if n == y else -(self.spaced_qmers[y].get_weight() - 1)
                        num_one_before_shift = self.spaced_qmers[y].get_weight() - 1
                        for init in range(starting_init, num_one_before_shift + 1):
                            temp = PreviousShiftMulti(Position([]), Position([]), i, init, n, 0)
                            for j in range(len(self.v_pos_ones[n])):
                                if (j - init) < len(self.v_pos_ones[y]) and (j - init) >= 0:
                                    if self.v_pos_ones[n][j] < i or self.v_pos_ones[y][j - init] != self.v_pos_ones[n][j] - i or self.v_pos_ones[y][j - init] not in pos_not_covered_yet:
                                        temp.one_to_remove.positions.append(j - init)
                                        temp.one_to_change.positions.append(j - init)
                                    else:
                                        temp.one_to_keep.positions.append(j - init)
                            if len(temp.one_to_keep.positions) > len(curr_best.one_to_keep.positions) or (len(temp.one_to_keep.positions) == len(curr_best.one_to_keep.positions) and temp.shift_min < curr_best.shift_min):
                                curr_best = temp

                for j in curr_best.one_to_keep.positions:
                    pos_not_covered_yet.remove(j)

                if len(curr_best.one_to_keep.positions) != 0:
                    self.multi_seed_info_col[y].group_previous[index].prev.append(curr_best)
                else:
                    done = True

            self.multi_seed_info_col[y].group_previous[index].not_covered = Position(pos_not_covered_yet)

        for j in range(len(self.multi_seed_info_col)):
            for k in range(len(self.multi_seed_info_col[j].group_previous[index].prev)):
                mask = 0
                for pos in self.multi_seed_info_col[j].group_previous[index].prev[k].one_to_keep.positions:
                    mask |= 3 << (pos * 2)
                self.multi_seed_info_col[j].group_previous[index].prev[k].mask = mask

    def set_multi_seed_info_row(self):
        self.multi_seed_info_row = MultiSeedInfoRow(info=[SeedInfo(group_previous=[groupPrevious(prev=[], not_covered=Position([]))], pos_ones=Position([])) for _ in range(self.get_length())], transient1=len(self.spaced_qs[0]) - 1, transient2=len(self.spaced_qs[0]) - 1)
        self.process_multi_seed_row(self.multi_seed_info_row.info, 0)

        furthest_pos_back = 0
        furthest_pos_front = 0
        for j in range(self.get_length()):
            for prev in self.multi_seed_info_row.info[j].group_previous[0].prev:
                furthest_pos_back = max(furthest_pos_back, prev.shift_min)
                furthest_pos_front = min(furthest_pos_front, prev.shift_min)

            self.multi_seed_info_row.info[j].pos_ones = Position(self.v_pos_ones[j])

        self.multi_seed_info_row.transient1 = furthest_pos_back
        self.multi_seed_info_row.transient2 = -furthest_pos_front

        for k in range(1, furthest_pos_back + 1):
            self.process_multi_seed_row(self.multi_seed_info_row.info, k)

        for k in range(furthest_pos_front, 0):
            self.process_multi_seed_row(self.multi_seed_info_row.info, k)

    def process_multi_seed_row(self, multi_seed_info_row, step):
        index = 0
        for y in range(self.get_length()):
            current_shift_group = groupPrevious(prev=[], not_covered=Position([]))
            current_shift_prev = []
            current_shift_group.prev = current_shift_prev
            multi_seed_info_row[y].group_previous.append(current_shift_group)
            index = len(multi_seed_info_row[y].group_previous) - 1
            pos_not_covered_yet = self.v_pos_pos_ones[y]

            done = False
            while not done:
                curr_best = PreviousShiftMulti(Position([]), Position([]), float('inf'))

                max_num_shifts = step - 1 if step > 0 else self.multi_seed_info_row.transient1
                min_num_shifts = -self.multi_seed_info_row.transient2 if step >= 0 else step + 1

                for n in range(y + 1):
                    starting_point = min_num_shifts if n < y else 1

                    for i in range(starting_point, max_num_shifts + 1):
                        if i >= 0:
                            starting_init = 1 if n == y else -(self.spaced_qmers[y].get_weight() - 1)
                            num_one_before_shift = self.spaced_qmers[y].get_weight() - 1
                            for init in range(starting_init, num_one_before_shift + 1):
                                temp = PreviousShiftMulti(Position([]), Position([]), i, init, n, 0)
                                for j in range(len(self.v_pos_ones[n])):
                                    if (j - init) < len(self.v_pos_ones[y]) and (j - init) >= 0:
                                        if self.v_pos_ones[n][j] < i or self.v_pos_ones[y][j - init] != self.v_pos_ones[n][j] - i or self.v_pos_ones[y][j - init] not in pos_not_covered_yet:
                                            temp.one_to_remove.positions.append(j - init)
                                            temp.one_to_change.positions.append(j - init)
                                        else:
                                            temp.one_to_keep.positions.append(j - init)
                                if len(temp.one_to_keep.positions) > len(curr_best.one_to_keep.positions) or (len(temp.one_to_keep.positions) == len(curr_best.one_to_keep.positions) and abs(temp.shift_min) < abs(curr_best.shift_min)):
                                    curr_best = temp
                        else:
                            num_one_after_shift = 0
                            while self.v_pos_ones[y][num_one_after_shift] < -i:
                                num_one_after_shift += 1

                            num_one_before_shift = len(self.v_pos_ones[y]) - 1
                            while self.v_pos_ones[y][num_one_before_shift] > -i:
                                num_one_before_shift -= 1
                            num_one_before_shift = len(self.v_pos_ones[y]) - num_one_before_shift
                            num_one_after_shift = -num_one_after_shift

                            starting_init = 1 if n == y else num_one_after_shift
                            for init in range(starting_init, num_one_before_shift + 1):
                                temp = PreviousShiftMulti(Position([]), Position([]), i, init, n, 0)
                                for j in range(len(self.v_pos_ones[n])):
                                    if (j - init) < len(self.v_pos_ones[y]) and (j - init) >= 0:
                                        if self.v_pos_ones[y][j - init] != self.v_pos_ones[n][j] - i or self.v_pos_ones[y][j - init] not in pos_not_covered_yet:
                                            temp.one_to_remove.positions.append(j - init)
                                            temp.one_to_change.positions.append(j - init)
                                        else:
                                            temp.one_to_keep.positions.append(j - init)
                                if len(temp.one_to_keep.positions) > len(curr_best.one_to_keep.positions) or (len(temp.one_to_keep.positions) == len(curr_best.one_to_keep.positions) and abs(temp.shift_min) < abs(curr_best.shift_min)):
                                    curr_best = temp

                for j in curr_best.one_to_keep.positions:
                    pos_not_covered_yet.remove(j)

                if len(curr_best.one_to_keep.positions) != 0:
                    multi_seed_info_row[y].group_previous[index].prev.append(curr_best)
                else:
                    done = True

            multi_seed_info_row[y].group_previous[index].not_covered = Position(pos_not_covered_yet)

        for j in range(len(multi_seed_info_row)):
            for k in range(len(multi_seed_info_row[j].group_previous[index].prev)):
                mask = 0
                for pos in multi_seed_info_row[j].group_previous[index].prev[k].one_to_keep.positions:
                    mask |= 3 << (pos * 2)
                multi_seed_info_row[j].group_previous[index].prev[k].mask = mask

def print_shift_multi(s: PreviousShiftMulti):
    print(f"one_to_remove= {s.one_to_remove.positions}")
    print(f"one_to_keep= {s.one_to_keep.positions}")
    print(f"one_exit= {s.one_exit}")
    print(f"shift_min= {s.shift_min}")
    print(f"seed_num= {s.seed_num}")
    print(f"mask= {bin(s.mask)}")
