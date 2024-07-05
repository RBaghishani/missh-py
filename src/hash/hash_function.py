from typing import List, Callable

# Define the necessary types
hash_type = int

class Hash_Err:
    def __init__(self):
        self.hash = 0
        self.errors = []

    def reset(self):
        self.hash = 0
        self.errors = []

    def push_back_error(self, pos):
        self.errors.append(pos)

    def sub_pos_err(self, shift, prev_hash):
        self.errors = [e - shift for e in prev_hash.errors if e - shift >= 0]

    def is_correct(self):
        return len(self.errors) == 0

    def size_error(self):
        return len(self.errors)

    def __getitem__(self, idx):
        return self.errors[idx]

    def sort_uniq_err(self):
        self.errors = sorted(set(self.errors))


Hash_Err_V = List[Hash_Err]

# Conversion functions
def char_to_int(ch: str) -> hash_type:
    return {'A': 0, 'C': 1, 'G': 2, 'T': 3}.get(ch, 4)  # 4 as error code

def char_to_int_complement(ch: str) -> hash_type:
    return {'A': 3, 'C': 2, 'G': 1, 'T': 0}.get(ch, 4)  # 4 as error code

# Print function
def print_hashes(v_hash_err: Hash_Err_V):
    for hash_err in v_hash_err:
        print(hash_err.hash)

# Hash functions
def get_hash(s_Str: str, startQmer: int, length: int, hash_err: Hash_Err, f_conversion: Callable[[str], hash_type]):
    hash_err.reset()
    for i in range(startQmer, startQmer + length):
        ch = f_conversion(s_Str[i])
        if ch == 4:
            hash_err.push_back_error(i)
        else:
            hash_err.hash |= ch << ((i - startQmer) * 2)

def get_hash_spaced(s_Str: str, startQmer: int, spaced_qmer: 'SpacedQmer', hash_err: Hash_Err, f_conversion: Callable[[str], hash_type]):
    hash_err.reset()
    pos_one = spaced_qmer.get_pos_one()
    for j in range(len(pos_one)):
        ch = f_conversion(s_Str[startQmer + pos_one[j]])
        if ch == 4:
            hash_err.push_back_error(j)
        else:
            hash_err.hash |= ch << (j * 2)

def get_hash_from_pos_one(s_Str: str, startQmer: int, pos_one: List[int], hash_err: Hash_Err, f_conversion: Callable[[str], hash_type]):
    hash_err.reset()
    for j in range(len(pos_one)):
        ch = f_conversion(s_Str[startQmer + pos_one[j]])
        if ch == 4:
            hash_err.push_back_error(j)
        else:
            hash_err.hash |= ch << (j * 2)

def get_hashes_speedup_previous(s_Str: str, length: int, vHash: Hash_Err_V, f_conversion: Callable[[str], hash_type]):
    vHash.clear()
    if len(s_Str) >= length:
        n_hashes = len(s_Str) - length + 1
        vHash.extend([Hash_Err() for _ in range(n_hashes)])

        get_hash(s_Str, 0, length, vHash[0], f_conversion)
        for pos in range(1, len(vHash)):
            prev_hash = vHash[pos - 1]
            curr_hash = vHash[pos]
            curr_hash.hash = prev_hash.hash >> 2
            curr_hash.sub_pos_err(1, prev_hash)

            enter = f_conversion(s_Str[pos + length - 1])
            if enter == 4:
                curr_hash.push_back_error(length - 1)
            else:
                curr_hash.hash |= enter << ((length - 1) * 2)

def get_hashes_naive(s_Str: str, spaced_qmer: 'SpacedQmer', vHash: Hash_Err_V, f_conversion: Callable[[str], hash_type]):
    vHash.clear()
    if len(s_Str) >= spaced_qmer.get_q():
        n_hashes = len(s_Str) - spaced_qmer.get_q() + 1
        vHash.extend([Hash_Err() for _ in range(n_hashes)])
        for pos in range(len(vHash)):
            get_hash_spaced(s_Str, pos, spaced_qmer, vHash[pos], f_conversion)

def compute_hash_for_speedup_previous(s_Str: str, pos_one_current: List[int], pos_one_prev: List[int], curr_sp_shift: 'PreviousShift',
                                      prev_hash_err: Hash_Err, idx_curr_hash: int, curr_hash_err: Hash_Err, f_conversion: Callable[[str], hash_type]):
    curr_hash_err.hash = prev_hash_err.hash >> (2 * curr_sp_shift.one_exit)
    if curr_sp_shift.one_to_remove:
        reset_one = 0
        for j in curr_sp_shift.one_to_remove:
            reset_one |= 3 << (j * 2)
        curr_hash_err.hash &= ~reset_one

    if not prev_hash_err.is_correct():
        for e in prev_hash_err.errors:
            curr_pos_one = e - curr_sp_shift.one_exit
            if curr_pos_one >= 0 and pos_one_prev[e] - curr_sp_shift.shift_min == pos_one_current[curr_pos_one]:
                curr_hash_err.push_back_error(curr_pos_one)

    for j in curr_sp_shift.one_to_change:
        i_to_change = j
        index_char = idx_curr_hash + pos_one_current[i_to_change]
        ch = f_conversion(s_Str[index_char])
        if ch == 4:
            curr_hash_err.push_back_error(i_to_change)
        else:
            curr_hash_err.hash |= ch << (i_to_change * 2)

    if len(pos_one_current) == len(pos_one_prev):
        for j in range(len(pos_one_current) - curr_sp_shift.one_exit, len(pos_one_current)):
            index_char = idx_curr_hash + pos_one_current[j]
            ch = f_conversion(s_Str[index_char])
            if ch == 4:
                curr_hash_err.push_back_error(j)
            else:
                curr_hash_err.hash |= ch << (j * 2)

    if not curr_hash_err.is_correct():
        curr_hash_err.sort_uniq_err()

def get_hashes_speedup_previous_spaced(s_Str: str, spaced_qmer: 'SpacedQmer', vHash: Hash_Err_V, f_conversion: Callable[[str], hash_type]):
    n_hashes = len(s_Str) - spaced_qmer.get_q() + 1
    vHash.clear()
    if n_hashes > 0:
        shift = spaced_qmer.get_shift_min_change()
        vHash.extend([Hash_Err() for _ in range(n_hashes)])
        get_hash_spaced(s_Str, 0, spaced_qmer, vHash[0], f_conversion)

        def get_hash(curr_idx_hash: int, curr_shift: 'PreviousShift'):
            curr_hash = vHash[curr_idx_hash]
            if spaced_qmer.get_weight() < curr_shift.get_size():
                get_hash_spaced(s_Str, curr_idx_hash, spaced_qmer, curr_hash, f_conversion)
            else:
                pos_hash_get = curr_idx_hash - curr_shift.shift_min
                prev_hash = vHash[pos_hash_get]
                compute_hash_for_speedup_previous(s_Str, spaced_qmer.get_pos_one(), spaced_qmer.get_pos_one(), curr_shift, prev_hash, curr_idx_hash, curr_hash, f_conversion)

        lim_max = len(vHash)
        lim_min = min(len(shift), lim_max)
        for i in range(1, lim_min):
            get_hash(i, shift[i])
        for i in range(lim_min, lim_max):
            get_hash(i, shift[-1])

def compute_hash_with_issh(curr_idx_hash: int, shifts: List['PreviousShift'], vHash: Hash_Err_V, pos_one: List[int], s_Str: str, f_conversion: Callable[[str], hash_type]):
    curr_hash = vHash[curr_idx_hash]
    curr_hash.hash = 0
    pos_not_covered_yet = shifts[0].one_to_change
    for shift in shifts:
        pos_hash_get = curr_idx_hash - shift.shift_min
        prev_hash = vHash[pos_hash_get]
        partial_hash = prev_hash.hash >> (2 * shift.one_exit)
        partial_hash &= shift.mask

        if not prev_hash.is_correct():
            for e in prev_hash.errors:
                curr_pos_one = e - shift.one_exit
                if curr_pos_one >= 0 and pos_one[e] - shift.shift_min == pos_one[curr_pos_one]:
                    curr_hash.push_back_error(curr_pos_one)

        curr_hash.hash |= partial_hash

    for i_to_change in pos_not_covered_yet:
        index_char = curr_idx_hash + pos_one[i_to_change]
        ch = f_conversion(s_Str[index_char])
        if ch == 4:
            curr_hash.push_back_error(i_to_change)
        else:
            curr_hash.hash |= ch << (i_to_change * 2)

def get_hashes_with_issh(s_Str: str, spaced_qmer: 'SpacedQmer', vHash: Hash_Err_V, f_conversion: Callable[[str], hash_type]):
    n_hashes = len(s_Str) - spaced_qmer.get_q() + 1
    vHash.clear()
    if n_hashes > 0:
        v_shifts = spaced_qmer.get_multiple_shifts()
        vHash.extend([Hash_Err() for _ in range(n_hashes)])
        get_hash_spaced(s_Str, 0, spaced_qmer, vHash[0], f_conversion)

        pos_one = spaced_qmer.get_pos_one()
        for i in range(1, min(n_hashes, len(v_shifts))):
            if not v_shifts[i]:
                get_hash_spaced(s_Str, i, spaced_qmer, vHash[i], f_conversion)
            else:
                compute_hash_with_issh(i, v_shifts[i], vHash, pos_one, s_Str, f_conversion)

        for i in range(len(v_shifts), n_hashes):
            compute_hash_with_issh(i, v_shifts[0], vHash, pos_one, s_Str, f_conversion)
