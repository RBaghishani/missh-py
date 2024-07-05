class SpacedQmer:
    def __init__(self, spaced_qmer: str = "", num_prev: int = 0):
        self.num_prev = num_prev
        self.spaced_q = spaced_qmer
        self.pos_one = self._get_pos_one(spaced_qmer)
        self.pos_pos_one = []
        self.shift_min_change = []
        self.multiple_shifts = []
        if spaced_qmer:
            self.reset(spaced_qmer, num_prev)

    def _get_pos_one(self, spaced_qmer: str):
        return [i for i, char in enumerate(spaced_qmer) if char == '1']

    def get_pos_one(self):
        return self.pos_one

    def get_q(self):
        return len(self.spaced_q)

    def get_weight(self):
        return len(self.pos_one)

    def reset(self, spaced_qmer, numprev):
        self.num_prev = numprev
        self.spaced_q = spaced_qmer
        self.save_index_one()
        self.shift_min_change = []
        self.get_shift_max(self.shift_min_change)
        self.set_all_multiple_shift()

    def save_index_one(self):
        self.pos_one = []
        self.pos_pos_one = []
        k = 0
        for i, char in enumerate(self.spaced_q):
            if self.is_one(i):
                self.pos_one.append(i)
                self.pos_pos_one.append(k)
                k += 1

    def is_one(self, index):
        return self.spaced_q[index] == '1'

    def get_shift_max(self, shift_max):
        shift_max.clear()
        shift_max.extend([{'one_to_remove': [], 'one_to_change': [], 'one_to_keep': [], 'one_exit': 0, 'shift_min': 0, 'mask': 0} for _ in range(len(self.spaced_q))])
        init = 0
        for i in range(1, len(self.spaced_q)):
            find = False
            for j in range(init, len(self.pos_one)):
                if self.pos_one[j] >= i:
                    init = j
                    find = True
                    break
            if not find:
                init = len(self.pos_one)
            for j in range(init, len(self.pos_one)):
                if self.pos_one[j - init] != self.pos_one[j] - i:
                    shift_max[i]['one_to_remove'].append(j - init)
                    shift_max[i]['one_to_change'].append(j - init)
                else:
                    shift_max[i]['one_to_keep'].append(j - init)
            shift_max[i]['one_exit'] = init
            shift_max[i]['shift_min'] = i
            if i > 1 and len(shift_max[shift_max[i - 1]['shift_min']]['one_to_keep']) < len(shift_max[i]['one_to_keep']):
                shift_max[i] = shift_max[i - 1]

    def set_all_multiple_shift(self):
        self.multiple_shifts = [[]]
        self.set_multiple_shifts(0)
        furthest_pos = self.multiple_shifts[0][0]['shift_min']
        for shift in self.multiple_shifts[0]:
            if furthest_pos < shift['shift_min']:
                furthest_pos = shift['shift_min']
        furthest_pos += 1
        self.multiple_shifts.extend([[] for _ in range(furthest_pos)])
        for k in range(1, furthest_pos):
            self.set_multiple_shifts(k)
        self.set_bit_masks()

    def set_multiple_shifts(self, index):
        self.multiple_shifts[index] = []
        pos_not_covered_yet = self.pos_pos_one[:]
        num_max_previous = (index if index > len(self.spaced_q) or index == 0 else len(self.spaced_q)) if self.num_prev == 0 else (self.num_prev if index > len(self.spaced_q) or index == 0 else min(index, self.num_prev))

        for k in range(num_max_previous):
            if len(pos_not_covered_yet) <= 1:
                break
            curr_best = {'one_to_remove': [], 'one_to_change': [], 'one_to_keep': [], 'one_exit': 0, 'shift_min': float('inf'), 'mask': 0}
            max_num_shifts = index if index != 0 else len(self.spaced_q)
            for i in range(1, max_num_shifts + 1):
                num_one_before_shift = 0
                while num_one_before_shift < len(self.pos_one) and self.pos_one[num_one_before_shift] < i:
                    num_one_before_shift += 1
                for init in range(1, num_one_before_shift + 1):
                    temp = {'one_to_remove': [], 'one_to_change': [], 'one_to_keep': [], 'one_exit': init, 'shift_min': i, 'mask': 0}
                    for j in range(init, len(self.pos_one)):
                        if (j - init) < len(self.pos_one):
                            if self.pos_one[j] < i or self.pos_one[j - init] != self.pos_one[j] - i or not self.is_contained(pos_not_covered_yet, self.pos_one, self.pos_one[j - init]):
                                temp['one_to_remove'].append(j - init)
                                temp['one_to_change'].append(j - init)
                            else:
                                temp['one_to_keep'].append(j - init)
                    if len(temp['one_to_keep']) > len(curr_best['one_to_keep']):
                        curr_best = temp

            for j in curr_best['one_to_keep']:
                self.delete_element(pos_not_covered_yet, j)
            if len(curr_best['one_to_keep']) != 0:
                self.multiple_shifts[index].append(curr_best)

        if len(self.multiple_shifts[index]) > 0:
            self.multiple_shifts[index][0]['one_to_change'] = pos_not_covered_yet

    def set_bit_masks(self):
        for shifts in self.multiple_shifts:
            for shift in shifts:
                for pos in shift['one_to_remove']:
                    shift['mask'] |= (3 << (pos * 2))
                shift['mask'] = ~shift['mask']

    def delete_element(self, pos, index):
        if index in pos:
            pos.remove(index)

    def is_contained(self, pointer, pos, index):
        return index in [pos[p] for p in pointer]

def print_shift(shift):
    print("one_to_change= ", shift['one_to_change'])
    print("one_to_remove= ", shift['one_to_remove'])
    print("one_to_keep= ", shift['one_to_keep'])
    print("one_exit= ", shift['one_exit'])
    print("shift_min= ", shift['shift_min'])
    print("mask= ", format(shift['mask'], '042b'))
    print()

def printp(pos):
    print(" ".join(map(str, pos)))
