# multi_spaced_qmer.py

from collections import namedtuple
import sys
import numpy as np

# Define necessary types
PreviousShiftMulti = namedtuple('PreviousShiftMulti', [
    'one_to_remove', 'one_to_keep', 'one_exit', 'shift_min', 'seed_num', 'mask'
])

class MultiSpacedQmer:
    def __init__(self, spaced_qmers):
        self.spaced_qmers = spaced_qmers
        self.spaced_qs = [''] * self.get_length()
        self.v_pos_ones = [[] for _ in range(self.get_length())]
        self.v_pos_pos_ones = [[] for _ in range(self.get_length())]
        self.multi_seed_info_col = []
        self.multi_seed_info_row = namedtuple('MultiSeedInfoRow', ['info', 'transient1', 'transient2'])

        length = spaced_qmers[0].get_q()
        weight = spaced_qmers[0].get_weight()

        for i in range(self.get_length()):
            if length != spaced_qmers[i].get_q() or weight != spaced_qmers[i].get_weight():
                print("\nThe seed in the group does not have the same length or weight\n", file=sys.stderr)
                sys.exit(1)

            self.spaced_qs[i] = spaced_qmers[i].to_string()
            k = 0
            for j in range(len(self.spaced_qs[i])):
                if spaced_qmers[i].is_one(j):
                    self.v_pos_ones[i].append(j)
                    self.v_pos_pos_ones[i].append(k)
                    k += 1

        self.set_multi_seed_info_row()
        self.v_pos_pos_ones = [[] for _ in range(self.get_length())]

        for i in range(self.get_length()):
            k = 0
            for j in range(len(self.spaced_qs[i])):
                if spaced_qmers[i].is_one(j):
                    self.v_pos_pos_ones[i].append(k)
                    k += 1

        self.set_multi_seed_info_col()

    def get_length(self):
        return len(self.spaced_qmers)

    def get_multi_seed_info_col(self):
        return self.multi_seed_info_col

    def get_multi_seed_info_row(self):
        return self.multi_seed_info_row

    def reset(self, spaced_qmers):
        self.__init__(spaced_qmers)

    def set_multi_seed_info_col(self):
        self.multi_seed_info_col = [{} for _ in range(self.get_length())]

        for j in range(self.get_length()):
            self.multi_seed_info_col[j]['group_previous'] = [{}]
            self.multi_seed_info_col[j]['pos_ones'] = self.v_pos_ones[j]

        self.process_multi_seed_col(0)

        furthest_pos = 0
        for j in range(self.get_length()):
            for i in range(len(self.multi_seed_info_col[j]['group_previous'][0]['prev'])):
                furthest_pos = max(furthest_pos, self.multi_seed_info_col[j]['group_previous'][0]['prev'][i]['shift_min'])

        furthest_pos += 1

        for j in range(self.get_length()):
            self.multi_seed_info_col[j]['group_previous'].extend([{}] * (furthest_pos - 1))

        for k in range(1, furthest_pos):
            self.process_multi_seed_col(k)

    def process_multi_seed_col(self, index):
        for y in range(self.get_length()):
            current_shift_group = []
            self.multi_seed_info_col[y]['group_previous'][index]['prev'] = current_shift_group
            pos_not_covered_yet = self.v_pos_pos_ones[y][:]

            done = False
            while not done:
                curr_best = PreviousShiftMulti([], [], 0, 0, 0, 0)

                max_num_shifts = index - 1 if index != 0 else len(self.spaced_qs[index]) - 1

                for n in range(self.get_length()):
                    starting_point = 0 if n < y else 1

                    for i in range(starting_point, max_num_shifts + 1):
                        starting_init = 1 if n == y else -(self.spaced_qmers[y].get_weight() - 1)
                        num_one_before_shift = self.spaced_qmers[y].get_weight() - 1

                        for init in range(starting_init, num_one_before_shift + 1):
                            temp = PreviousShiftMulti([], [], 0, 0, 0, 0)

                            for j in range(len(self.v_pos_ones[n])):
                                if (j - init) < len(self.v_pos_ones[y]) and (j - init) >= 0:
                                    if self.v_pos_ones[n][j] < i or self.v_pos_ones[y][j - init] != self.v_pos_ones[n][j] - i or not self.is_contained(pos_not_covered_yet, self.v_pos_ones[y], self.v_pos_ones[y][j - init]):
                                        temp.one_to_remove.append(j - init)
                                    else:
                                        temp.one_to_keep.append(j - init)

                            temp = temp._replace(one_exit=init, shift_min=i, seed_num=n)
                            if len(temp.one_to_keep) > len(curr_best.one_to_keep) or (len(temp.one_to_keep) == len(curr_best.one_to_keep) and temp.shift_min < curr_best.shift_min):
                                curr_best = temp

                for j in curr_best.one_to_keep:
                    pos_not_covered_yet.remove(j)

                if len(curr_best.one_to_keep) != 0:
                    self.multi_seed_info_col[y]['group_previous'][index]['prev'].append(curr_best)
                else:
                    done = True

            self.multi_seed_info_col[y]['group_previous'][index]['not_covered'] = pos_not_covered_yet

        for j in range(len(self.multi_seed_info_col)):
            for k in range(len(self.multi_seed_info_col[j]['group_previous'][index]['prev'])):
                for i in range(len(self.multi_seed_info_col[j]['group_previous'][index]['prev'][k].one_to_keep)):
                    self.multi_seed_info_col[j]['group_previous'][index]['prev'][k].mask |= 3 << (self.multi_seed_info_col[j]['group_previous'][index]['prev'][k].one_to_keep[i] * 2)

    def set_multi_seed_info_row(self):
        self.multi_seed_info_row.info = [{} for _ in range(self.get_length())]
        self.multi_seed_info_row.transient1 = len(self.spaced_qs[0]) - 1
        self.multi_seed_info_row.transient2 = len(self.spaced_qs[0]) - 1

        self.process_multi_seed_row(self.multi_seed_info_row.info, 0)

        furthest_pos_back = 0
        furthest_pos_front = 0

        for j in range(self.get_length()):
            for i in range(len(self.multi_seed_info_row.info[j]['group_previous'][0]['prev'])):
                furthest_pos_back = max(furthest_pos_back, self.multi_seed_info_row.info[j]['group_previous'][0]['prev'][i]['shift_min'])
                furthest_pos_front = min(furthest_pos_front, self.multi_seed_info_row.info[j]['group_previous'][0]['prev'][i]['shift_min'])

            self.multi_seed_info_row.info[j]['pos_ones'] = self.v_pos_ones[j]

        self.multi_seed_info_row.transient1 = furthest_pos_back
        self.multi_seed_info_row.transient2 = -furthest_pos_front

        for k in range(1, furthest_pos_back + 1):
            self.process_multi_seed_row(self.multi_seed_info_row.info, k)

        for k in range(furthest_pos_front, 0):
            self.process_multi_seed_row(self.multi_seed_info_row.info, k)

    def process_multi_seed_row(self, multi_seed_info_row, step):
        index = 0
        for y in range(self.get_length()):
            current_shift_group = {}
            current_shift_prev = []
            current_shift_group['prev'] = current_shift_prev
            multi_seed_info_row[y]['group_previous'].append(current_shift_group)
            index = len(multi_seed_info_row[y]['group_previous']) - 1
            pos_not_covered_yet = self.v_pos_pos_ones[y][:]

            done = False
            while not done:
                curr_best = PreviousShiftMulti([], [], 0, 0, 0, 0)
                max_num_shifts = step - 1 if step > 0 else self.multi_seed_info_row.transient1
                min_num_shifts = -self.multi_seed_info_row.transient2 if step >= 0 else step + 1

                for n in range(y + 1):
                    starting_point = min_num_shifts if n < y else 1

                    for i in range(starting_point, max_num_shifts + 1):
                        if i >= 0:
                            starting_init = 1 if n == y else -(self.spaced_qmers[y].get_weight() - 1)
                            num_one_before_shift = self.spaced_qmers[y].get_weight() - 1

                            for init in range(starting_init, num_one_before_shift + 1):
                                temp = PreviousShiftMulti([], [], 0, 0, 0, 0)

                                for j in range(len(self.v_pos_ones[n])):
                                    if (j - init) < len(self.v_pos_ones[y]) and (j - init) >= 0:
                                        if self.v_pos_ones[n][j] < i or self.v_pos_ones[y][j - init] != self.v_pos_ones[n][j] - i or not self.is_contained(pos_not_covered_yet, self.v_pos_ones[y], self.v_pos_ones[y][j - init]):
                                            temp.one_to_remove.append(j - init)
                                        else:
                                            temp.one_to_keep.append(j - init)

                                temp = temp._replace(one_exit=init, shift_min=i, seed_num=n)
                                if len(temp.one_to_keep) > len(curr_best.one_to_keep) or (len(temp.one_to_keep) == len(curr_best.one_to_keep) and abs(temp.shift_min) < abs(curr_best.shift_min)):
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
                                temp = PreviousShiftMulti([], [], 0, 0, 0, 0)

                                for j in range(len(self.v_pos_ones[n])):
                                    if (j - init) < len(self.v_pos_ones[y]) and (j - init) >= 0:
                                        if self.v_pos_ones[y][j - init] != self.v_pos_ones[n][j] - i or not self.is_contained(pos_not_covered_yet, self.v_pos_ones[y], self.v_pos_ones[y][j - init]):
                                            temp.one_to_remove.append(j - init)
                                        else:
                                            temp.one_to_keep.append(j - init)

                                temp = temp._replace(one_exit=init, shift_min=i, seed_num=n)
                                if len(temp.one_to_keep) > len(curr_best.one_to_keep) or (len(temp.one_to_keep) == len(curr_best.one_to_keep) and abs(temp.shift_min) < abs(curr_best.shift_min)):
                                    curr_best = temp

                for j in curr_best.one_to_keep:
                    pos_not_covered_yet.remove(j)

                if len(curr_best.one_to_keep) != 0:
                    multi_seed_info_row[y]['group_previous'][index]['prev'].append(curr_best)
                else:
                    done = True

            multi_seed_info_row[y]['group_previous'][index]['not_covered'] = pos_not_covered_yet

        for j in range(len(multi_seed_info_row)):
            for k in range(len(multi_seed_info_row[j]['group_previous'][index]['prev'])):
                for i in range(len(multi_seed_info_row[j]['group_previous'][index]['prev'][k].one_to_keep)):
                    multi_seed_info_row[j]['group_previous'][index]['prev'][k].mask |= 3 << (multi_seed_info_row[j]['group_previous'][index]['prev'][k].one_to_keep[i] * 2)

    def is_contained(self, pos_not_covered_yet, v_pos_ones, value):
        return value in pos_not_covered_yet

    def delete_element(self, pos_not_covered_yet, value):
        if value in pos_not_covered_yet:
            pos_not_covered_yet.remove(value)

    def print_shift_multi(self, s):
        print("\nC_", s.shift_min - s.one_exit, ",", s.shift_min)
        print("\none_to_remove= ", s.one_to_remove)
        print("\none_to_keep= ", s.one_to_keep)
        print("\none_exit= ", s.one_exit)
        print("\nshift_min= ", s.shift_min)
        print("\nseed_num= ", s.seed_num)
        if s.mask != 0:
            print("\nMask: ", format(s.mask, '042b'))
        print("\n\n")

