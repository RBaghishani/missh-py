from typing import List, Callable
from hash_type import Hash_Err, hash_type, Hash_Err_V, Hash_Err_V_V
from spaced.spaced_qmer import SpacedQmer, SpacedQmer_Multi, PreviousShiftMulti, V_PreviusShiftMulti
from spaced.multi_spaced_qmer import MultiSeedInfo, MultiSeedInfoRow
import math
import numpy as np
from concurrent.futures import ThreadPoolExecutor

def get_hash(s_Str: str, startQmer: int, length: int, hash_err: Hash_Err, f_conversion: Callable[[str], int]):
    hash_err.reset()
    for i in range(startQmer, startQmer + length):
        ch = f_conversion(s_Str[i])
        if ch == 4:
            hash_err.push_back_error(i)
        else:
            hash_err.hash |= ch << ((i - startQmer) * 2)

def get_hash_spaced(s_Str: str, startQmer: int, spaced_qmer: SpacedQmer, hash_err: Hash_Err, f_conversion: Callable[[str], int]):
    hash_err.reset()
    pos_one = spaced_qmer.get_pos_one()
    for j in range(len(pos_one)):
        ch = f_conversion(s_Str[startQmer + pos_one[j]])
        if ch == 4:
            hash_err.push_back_error(j)
        else:
            hash_err.hash |= ch << (j * 2)

def get_hash_from_pos_one(s_Str: str, startQmer: int, pos_one: List[int], hash_err: Hash_Err, f_conversion: Callable[[str], int]):
    hash_err.reset()
    for j in range(len(pos_one)):
        ch = f_conversion(s_Str[startQmer + pos_one[j]])
        if ch == 4:
            hash_err.push_back_error(j)
        else:
            hash_err.hash |= ch << (j * 2)

def compute_hash_for_speedup_previous(s_Str: str, pos_one_current: List[int], pos_one_prev: List[int], curr_sp_shift: PreviousShiftMulti, prev_hash_err: Hash_Err, idx_curr_hash: int, curr_hash_err: Hash_Err, f_conversion: Callable[[str], int]):
    curr_hash_err.hash = prev_hash_err.hash
    curr_hash_err.hash >>= 2 * curr_sp_shift.one_exit
    if curr_sp_shift.one_to_remove:
        reset_one = 0
        for j in curr_sp_shift.one_to_remove:
            reset_one |= (hash_type) << (j * 2)
        curr_hash_err.hash &= ~reset_one
    if not prev_hash_err.is_correct():
        for e in range(prev_hash_err.size_error()):
            curr_pos_one = prev_hash_err[e] - curr_sp_shift.one_exit
            if curr_pos_one >= 0 and pos_one_prev[prev_hash_err[e]] - curr_sp_shift.shift_min == pos_one_current[curr_pos_one]:
                curr_hash_err.push_back_error(curr_pos_one)
    for j in curr_sp_shift.one_to_change:
        i_to_change = curr_sp_shift.one_to_change[j]
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

def get_hashes_speedup_multi_previous(s_Str: str, spaced_qmers: SpacedQmer_Multi, vHashes: Hash_Err_V_V, f_conversion: Callable[[str], int]):
    def get_hash(curr_spaced, curr_idx_hash, curr_shift):
        curr_hash = vHashes[curr_spaced][curr_idx_hash]
        if curr_shift.current_sp_ptr.get_weight() < curr_shift.get_size():
            get_hash_spaced(s_Str, curr_idx_hash, curr_shift.current_sp_ptr, curr_hash, f_conversion)
        else:
            pos_hash_get = curr_idx_hash - curr_shift.shift_min
            prev_hash = vHashes[curr_shift.prev_qmer][pos_hash_get]
            compute_hash_for_speedup_previous(s_Str, curr_shift.current_sp_ptr.get_pos_one(), curr_shift.prev_sp_ptr.get_pos_one(), curr_shift, prev_hash, curr_idx_hash, curr_hash, f_conversion)
    n_hashes = len(s_Str) - spaced_qmers[0].get_q() + 1
    vHashes.clear()
    if n_hashes > 0:
        v_shift_min = spaced_qmers.get_shift_min()
        vHashes.extend([[Hash_Err() for _ in range(n_hashes)] for _ in range(spaced_qmers.size())])
        for s in range(spaced_qmers.size()):
            if not v_shift_min[s][0].is_correct_spaced_previous():
                get_hash_spaced(s_Str, 0, spaced_qmers[s], vHashes[s][0], f_conversion)
            else:
                get_hash(s, 0, v_shift_min[s][0])
        lim_max = len(vHashes[0])
        lim_min = min(len(v_shift_min[0]), lim_max)
        for i in range(1, lim_min):
            for ss in range(spaced_qmers.size()):
                get_hash(ss, i, v_shift_min[ss][i])
        for i in range(lim_min, lim_max):
            for ss in range(spaced_qmers.size()):
                get_hash(ss, i, v_shift_min[ss][-1])
    else:
        vHashes.extend([[] for _ in range(spaced_qmers.size())])

def get_hashes_with_issh_multi_v1(s_Str: str, spaced_qmers: List[SpacedQmer], VV_shifts: List[List[List[PreviousShiftMulti]]], v_pos_one: List[List[int]], max_transient_length: int, vvHash: Hash_Err_V_V, f_conversion: Callable[[str], int]):
    n_hashes = len(s_Str) - spaced_qmers[0].get_q() + 1
    if n_hashes > 0:
        for j in range(len(spaced_qmers)):
            vvHash[j].clear()
            vvHash[j].extend([Hash_Err() for _ in range(n_hashes)])
            curr_seed = VV_shifts[j]
            curr_max = len(curr_seed)
            get_hash_spaced(s_Str, 0, spaced_qmers[j], vvHash[j][0], f_conversion)
            i = 1
            while i < n_hashes and i < max_transient_length:
                if not curr_seed[i]:
                    get_hash_spaced(s_Str, i, spaced_qmers[j], vvHash[j][i], f_conversion)
                else:
                    shift_group = 0 if i >= curr_max else i
                    compute_hash_with_issh(i, curr_seed[shift_group], vvHash[j], v_pos_one[j], s_Str, f_conversion)
                i += 1
        for i in range(i, n_hashes):
            compute_hash_with_issh_multi_v1(i, VV_shifts, vvHash, v_pos_one, s_Str, f_conversion)
    else:
        for j in range(len(VV_shifts)):
            vvHash[j].clear()

def compute_hash_with_issh_multi_v1(curr_idx_hash: int, VV_shifts: List[List[List[PreviousShiftMulti]]], vvHash: Hash_Err_V_V, v_pos_one: List[List[int]], s_Str: str, f_conversion: Callable[[str], int]):
    pos_not_covered_yet = VV_shifts[0][0][0].one_to_change
    i_to_change = pos_not_covered_yet[0]
    index_char = curr_idx_hash + v_pos_one[0][i_to_change]
    ch = f_conversion(s_Str[index_char])
    for j in range(len(VV_shifts)):
        curr_hash = vvHash[j][curr_idx_hash]
        curr_hash.hash = 0
        for k in range(len(VV_shifts[j][0])):
            pos_hash_get = curr_idx_hash - VV_shifts[j][0][k].shift_min
            prev_hash = vvHash[j][pos_hash_get]
            partial_hash = prev_hash.hash
            curr_sp_shift = VV_shifts[j][0][k]
            partial_hash >>= 2 * curr_sp_shift.one_exit
            partial_hash &= curr_sp_shift.mask
            if not prev_hash.is_correct():
                pass
            curr_hash.hash |= partial_hash
        if ch == 4:
            curr_hash.push_back_error(i_to_change)
        else:
            curr_hash.hash |= ch << (i_to_change * 2)

def compute_hash_with_issh(curr_idx_hash: int, shifts: List[PreviousShiftMulti], vHash: Hash_Err_V, pos_one: List[int], s_Str: str, f_conversion: Callable[[str], int]):
    curr_hash = vHash[curr_idx_hash]
    curr_hash.hash = 0
    pos_not_covered_yet = shifts[0].one_to_change
    for k in range(len(shifts)):
        pos_hash_get = curr_idx_hash - shifts[k].shift_min
        prev_hash = vHash[pos_hash_get]
        partial_hash = prev_hash.hash
        curr_sp_shift = shifts[k]
        partial_hash >>= 2 * curr_sp_shift.one_exit
        partial_hash &= curr_sp_shift.mask
        if not prev_hash.is_correct():
            pass
        curr_hash.hash |= partial_hash
    for i_to_change in pos_not_covered_yet:
        index_char = curr_idx_hash + pos_one[i_to_change]
        ch = f_conversion(s_Str[index_char])
        if ch == 4:
            curr_hash.push_back_error(i_to_change)
        else:
            curr_hash.hash |= ch << (i_to_change * 2)

def get_hashes_with_issh_multi_col(s_Str: str, VV_shifts: MultiSeedInfo, vvHash: Hash_Err_V_V, f_conversion: Callable[[str], int]):
    n_hashes = len(s_Str) - VV_shifts[0].pos_ones[-1]
    if n_hashes > 0:
        for j in range(len(VV_shifts)):
            vvHash[j].clear()
            vvHash[j].extend([Hash_Err() for _ in range(n_hashes)])
        for i in range(n_hashes):
            compute_hash_with_issh_multi_col(i, VV_shifts, vvHash, s_Str, f_conversion)
    else:
        for j in range(len(VV_shifts)):
            vvHash[j].clear()

def compute_hash_with_issh_multi_col(curr_idx_hash: int, VV_shifts: MultiSeedInfo, vvHash: Hash_Err_V_V, s_Str: str, f_conversion: Callable[[str], int]):
    shift_group = curr_idx_hash < len(VV_shifts[0].group_previous) - 1 if curr_idx_hash + 1 else 0
    for j in range(len(VV_shifts)):
        curr_hash = vvHash[j][curr_idx_hash]
        curr_hash.hash = 0
        curr_group = VV_shifts[j].group_previous[shift_group]
        curr_group_prev = curr_group.prev
        if not VV_shifts[j].group_previous[shift_group].prev:
            get_hash_from_pos_one(s_Str, curr_idx_hash, VV_shifts[j].pos_ones, vvHash[j][curr_idx_hash], f_conversion)
        else:
            for k in range(len(curr_group_prev)):
                curr_sp_shift = curr_group_prev[k]
                pos_hash_get = curr_idx_hash - curr_sp_shift.shift_min
                prev_hash = vvHash[curr_sp_shift.seed_num][pos_hash_get]
                partial_hash = prev_hash.hash
                if curr_sp_shift.one_exit >= 0:
                    partial_hash >>= 2 * curr_sp_shift.one_exit
                else:
                    partial_hash <<= -2 * curr_sp_shift.one_exit
                partial_hash &= curr_sp_shift.mask
                if not prev_hash.is_correct():
                    pass
                curr_hash.hash |= partial_hash
            pos_not_covered_yet = curr_group.not_covered
            pos_one = VV_shifts[j].pos_ones
            for i_to_change in pos_not_covered_yet:
                ch = f_conversion(s_Str[curr_idx_hash + pos_one[i_to_change]])
                if ch == 4:
                    curr_hash.push_back_error(i_to_change)
                else:
                    curr_hash.hash |= ch << (i_to_change * 2)

def get_hashes_with_issh_multi_col_parallel(s_Str: str, VV_shifts: MultiSeedInfo, vvHash: Hash_Err_V_V, f_conversion: Callable[[str], int]):
    n_hashes = len(s_Str) - VV_shifts[0].pos_ones[-1]
    if n_hashes > 0:
        for j in range(len(VV_shifts)):
            vvHash[j].clear()
            vvHash[j].extend([Hash_Err() for _ in range(n_hashes)])
        with ThreadPoolExecutor() as executor:
            for i in range(n_hashes):
                executor.submit(compute_hash_with_issh_multi_col_parallel, i, i, VV_shifts, vvHash, s_Str, f_conversion)
    else:
        for j in range(len(VV_shifts)):
            vvHash[j].clear()

def compute_hash_with_issh_multi_col_parallel(curr_idx_hash: int, offset: int, VV_shifts: MultiSeedInfo, vvHash: Hash_Err_V_V, s_Str: str, f_conversion: Callable[[str], int]):
    shift_group = offset < len(VV_shifts[0].group_previous) - 1 if offset + 1 else 0
    for j in range(len(VV_shifts)):
        curr_hash = vvHash[j][curr_idx_hash]
        curr_hash.hash = 0
        curr_group = VV_shifts[j].group_previous[shift_group]
        curr_group_prev = curr_group.prev
        if not VV_shifts[j].group_previous[shift_group].prev:
            get_hash_from_pos_one(s_Str, curr_idx_hash, VV_shifts[j].pos_ones, vvHash[j][curr_idx_hash], f_conversion)
        else:
            for k in range(len(curr_group_prev)):
                curr_sp_shift = curr_group_prev[k]
                pos_hash_get = curr_idx_hash - curr_sp_shift.shift_min
                prev_hash = vvHash[curr_sp_shift.seed_num][pos_hash_get]
                partial_hash = prev_hash.hash
                if curr_sp_shift.one_exit >= 0:
                    partial_hash >>= 2 * curr_sp_shift.one_exit
                else:
                    partial_hash <<= -2 * curr_sp_shift.one_exit
                partial_hash &= curr_sp_shift.mask
                if not prev_hash.is_correct():
                    pass
                curr_hash.hash |= partial_hash
            pos_not_covered_yet = curr_group.not_covered
            pos_one = VV_shifts[j].pos_ones
            for i_to_change in pos_not_covered_yet:
                ch = f_conversion(s_Str[curr_idx_hash + pos_one[i_to_change]])
                if ch == 4:
                    curr_hash.push_back_error(i_to_change)
                else:
                    curr_hash.hash |= ch << (i_to_change * 2)

def get_hashes_with_issh_multi_row(s_Str: str, rowInfo: MultiSeedInfoRow, vvHash: Hash_Err_V_V, f_conversion: Callable[[str], int]):
    VV_shifts = rowInfo.info
    n_hashes = len(s_Str) - VV_shifts[0].pos_ones[-1]
    transient1 = rowInfo.transient1
    transient2 = rowInfo.transient2
    if n_hashes > 0 and n_hashes >= (transient1 + transient2):
        for j in range(len(VV_shifts)):
            vvHash[j].clear()
            vvHash[j].extend([Hash_Err() for _ in range(n_hashes)])
            V_shift = VV_shifts[j].group_previous
            for curr_idx_hash in range(n_hashes):
                shift_group = (curr_idx_hash >= transient1 and curr_idx_hash < n_hashes - transient2) and 0 or (curr_idx_hash < transient1 and (curr_idx_hash + 1) or (len(V_shift) - (n_hashes - curr_idx_hash)))
                curr_group = V_shift[shift_group]
                curr_group_prev = curr_group.prev
                curr_hash = vvHash[j][curr_idx_hash]
                curr_hash.hash = 0
                if not VV_shifts[j].group_previous[shift_group].prev:
                    get_hash_from_pos_one(s_Str, curr_idx_hash, VV_shifts[j].pos_ones, vvHash[j][curr_idx_hash], f_conversion)
                else:
                    for k in range(len(curr_group_prev)):
                        curr_sp_shift = curr_group_prev[k]
                        pos_hash_get = curr_idx_hash - curr_sp_shift.shift_min
                        prev_hash = vvHash[curr_sp_shift.seed_num][pos_hash_get]
                        partial_hash = prev_hash.hash
                        if curr_sp_shift.one_exit >= 0:
                            partial_hash >>= 2 * curr_sp_shift.one_exit
                        else:
                            partial_hash <<= -2 * curr_sp_shift.one_exit
                        partial_hash &= curr_sp_shift.mask
                        if not prev_hash.is_correct():
                            pass
                        curr_hash.hash |= partial_hash
                    pos_not_covered_yet = curr_group.not_covered
                    pos_one = VV_shifts[j].pos_ones
                    for i_to_change in pos_not_covered_yet:
                        ch = f_conversion(s_Str[curr_idx_hash + pos_one[i_to_change]])
                        if ch == 4:
                            curr_hash.push_back_error(i_to_change)
                        else:
                            curr_hash.hash |= ch << (i_to_change * 2)
    elif n_hashes < (transient1 + transient2) and n_hashes > 0:
        for j in range(len(VV_shifts)):
            vvHash[j].clear()
            vvHash[j].extend([Hash_Err() for _ in range(n_hashes)])
            for curr_idx_hash in range(n_hashes):
                get_hash_from_pos_one(s_Str, curr_idx_hash, VV_shifts[j].pos_ones, vvHash[j][curr_idx_hash], f_conversion)
    else:
        for j in range(len(VV_shifts)):
            vvHash[j].clear()