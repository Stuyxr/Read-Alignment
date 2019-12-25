import random
from time import time


def is_lms(t, i):
    if i == 0:
        return t[i] == 'S'
    else:
        return t[i-1] == 'L' and t[i] == 'S'


def classify_ls(S):
    t = [None] * len(S)
    t[-1] = 'S'
    for i in reversed(range(len(S)-1)):
        if S[i] < S[i+1]:
            t[i] = 'S'
        elif S[i] > S[i+1]:
            t[i] = 'L'
        else:
            t[i] = t[i+1]
    lmss = [i for i in range(len(S)) if is_lms(t, i)]
    return t, lmss


def different_lms_substring(S, t, i, j):
    for k in range(len(S)):
        if S[i+k] != S[j+k]:
            return True
        elif (k > 0) and (is_lms(t, i+k) or is_lms(t, j+k)):
            return False
    raise RuntimeError


def induced_sort(S, k, t, lmss):
    sa = [None] * len(S)

    # 文字cはbins[c]からbins[c+1]-1の間に入る
    bins = [0] * k
    for c in S:
        bins[c] += 1
    for i in range(1, k):
        bins[i] = bins[i-1] + bins[i]
    bins.insert(0, 0)
    # print(bins)

    # step1 LMSをビンの後ろから詰める
    counts = [0] * k
    for i in lmss:
        c = S[i]
        sa[bins[c+1]-1-counts[c]] = i
        counts[c] += 1
        # print(sa)

    # step2 L型を詰める
    counts = [0] * k
    for i in sa:
        if i == None or i == 0 or t[i-1] == 'S':
            continue
        c = S[i-1]
        sa[bins[c]+counts[c]] = i-1
        counts[c] += 1
        # print(sa)

    # step3 S型を詰める
    counts = [0] * k
    for i in reversed(sa):
        if i == None or i == 0 or t[i-1] == 'L':
            continue
        c = S[i-1]
        sa[bins[c+1]-1-counts[c]] = i-1
        counts[c] += 1
        # print(sa)

    return sa


def sa_is(S, k):
    if type(S) is str:
        S = list(map(ord, S))
    t, lmss = classify_ls(S)

    # induced sort 1回目
    sa = induced_sort(S, k, t, lmss)

    # LMSのsuffixだけ取得
    sa = [i for i in sa if is_lms(t, i)]
    # print(sa)

    # LMSに番号を振る
    nums = [None] * len(S)
    nums[sa[0]] = 0
    num = 0
    for i, j in zip(sa, sa[1:]):
        if different_lms_substring(S, t, i, j):
            num += 1
        nums[j] = num
    nums = [n for n in nums if n is not None]

    # 再帰的に並べ替える
    if num < len(nums)-1:
        sa = sa_is(nums, num+1)
    else:
        sa = [None] * len(nums)
        for i, n in enumerate(nums):
            sa[n] = i

    # 最後に再度induced sort
    seed = [lmss[i] for i in sa]
    sa = induced_sort(S, k, t, seed)
    return sa
