import loading
from time import time
from random import randint

dict1 = {'$': 0, 'A': 1, 'C': 2, 'G': 3, 'T': 4, 'Z': 5}
dict2 = ['$', 'A', 'C', 'G', 'T', 'Z']
k = 256
# GGAACTCTCCCGCACCTTCGCTCACGTTGAT

def get_bwt_array(s, rev=0):
    """
    C[c]的定义是，字典序小于字符c的所有字符个数。
    Occ[r, c]表示在bwt[]中第r行之前出现字符c的个数。
    :param s: 参考基因序列s
    :param rev: 是否为反向计算
    :return: bwt数组，C数组，Occ数组，后缀数组
    """
    n = len(s)
    bwt = []
    C = [0] * k
    if rev == 0:
        with open('data/sa.txt', 'r') as f:
            suffix_arr = f.read()
            suffix_arr = [int(x) for x in suffix_arr.split(' ')]
    else:
        with open('data/rev_sa.txt', 'r') as f:
            suffix_arr = f.read()
            suffix_arr = [int(x) for x in suffix_arr.split(' ')]
    for i in range(len(s)):
        for x in dict1.keys():
            if s[i] < x:
                C[dict1[x]] += 1
    for i in range(n):
        bwt.append(s[suffix_arr[i] - 1])
    Occ = []
    for i in range(len(bwt)):
        if i == 0:
            Occ.append([0, 0, 0, 0, 0, 0])
        else:
            Occ.append(list(Occ[i - 1]))
        Occ[i][dict1[bwt[i]]] += 1
    Occ.append(list(Occ[-1]))
    for i in range(len(Occ)-1):
        Occ[i][dict1[bwt[i]]] -= 1
    return bwt, C, Occ, suffix_arr


def lfc(r, c):
    """
    lfc函数
    :param r: 行数r
    :param c: 字符c
    :return: lfc值
    """
    c = dict1[c]
    return C[c] + Occ[r][c]


def r_lfc(r, c):
    """
    反向lfc函数
    :param r: 行数r
    :param c: 字符c
    :return: 反向lfc值
    """
    c = dict1[c]
    return r_C[c] + r_Occ[r][c]


def exact_match(P):
    """
    字符串精确匹配
    :param P: 字符串P
    :return: 匹配上的集合
    """
    p = len(P) - 1
    sp = 0
    ep = len(ref)
    i = p
    while sp < ep and i >= 0:
        c = P[i]
        sp = lfc(sp, c)
        ep = lfc(ep, c)
        i = i - 1
    return set([x for x in range(sp, ep)])


def calculate_d(P):
    """
    计算错误上限
    :param P: read串
    :return: 数组d
    """
    ref_length = len(r_ref)
    sp = 0
    ep = ref_length
    z = 0
    d = [0] * len(P)
    for i in range(len(P)):
        c = P[i]
        sp = r_lfc(sp, c)
        ep = r_lfc(ep, c)
        if sp >= ep:
            z += 1
            sp = 0
            ep = ref_length
        d[i] = z
    return d


def inex_recur(P, i, z, sp, ep):
    """
    模糊匹配的搜索主体
    :param P: read串
    :param i: 当前搜索位置
    :param z: 错误上限
    :param sp: 左边界
    :param ep: 右边界
    :return: 符合条件的集合
    """
    # print(i, z, sp, ep, P[i])
    if i < 0:
        return set([x for x in range(sp, ep)])
    if z < d[i]:
        return set()
    s = set()
    # s = s | inex_recur(P, i-1, z-1, sp, ep)
    spp, epp = sp, ep
    for c in dict2[1:5]:
        sp = lfc(sp, c)
        ep = lfc(ep, c)
        if sp < ep:
            s = s | inex_recur(P, i, z-1, sp, ep)
            if c == P[i]:
                s = s | inex_recur(P, i-1, z, sp, ep)
            else:
                s = s | inex_recur(P, i-1, z-1, sp, ep)
        sp = spp
        ep = epp
    return s


def inexact_match(P, z):
    """
    模糊匹配
    :param P: read串
    :param z: 错误上限
    :return: 满足条件的集合
    """
    global d
    d = calculate_d(P)
    # print(d)
    return inex_recur(P, len(P)-1, z, 0, len(ref))


if __name__ == '__main__':
    ref = loading.load_ref('data/NC_008253.fna') + '$'
    # ref = ref[4910000:]
    # print(len(ref))
    with open('file.txt', 'w') as f:
        f.write(ref)
    reads1, reads2 = loading.load_reads('data/NC_008253_1.fastq', 'data/NC_008253_2.fastq')
    # st = time()
    bwt, C, Occ, sa = get_bwt_array(ref)
    # ed = time()
    # print(str(ed - st) + 'ms')
    # # print(exact_match('CATCAT'))
    # r_ref, r_C, r_Occ = str(), [], []
    #
    r_ref = ref[-2::-1] + '$'
    _, r_C, r_Occ, _ = get_bwt_array(r_ref)
    # # pos = inexact_match('CGG', 1)
    # #
    loading.write_head('result.sam')
    for i in range(1000):
        read1 = reads1[i]
        read2 = reads2[i]
        l = len(read1)
        if len(set(read1)) > 4 or len(set(read2)) > 4:
            continue
        read3 = str()
        read4 = str()
        for ch in read1:
            read3 += dict2[5-dict1[ch]]
        read3 = read3[::-1]
        for ch in read2:
            read4 += dict2[5-dict1[ch]]
        read4 = read4[::-1]
        # 以下位准确匹配
        pos_set = exact_match(read1)
        pos_list1 = [sa[x] for x in pos_set]
        pos_set = exact_match(read2)
        pos_list2 = [sa[x] for x in pos_set]
        pos_set = exact_match(read3)
        pos_list3 = [sa[x] for x in pos_set]
        pos_set = exact_match(read4)
        pos_list4 = [sa[x] for x in pos_set]
        # 以下位模糊匹配
        # pos_set = inexact_match(read1, 12)
        # pos_list1 = [sa[x] for x in pos_set]
        # pos_set = inexact_match(read2, 12)
        # pos_list2 = [sa[x] for x in pos_set]
        # pos_set = inexact_match(read3, 12)
        # pos_list3 = [sa[x] for x in pos_set]
        # pos_set = inexact_match(read4, 12)
        # pos_list4 = [sa[x] for x in pos_set]
        loading.write_sam('result.sam', '@NC_008253.fastq.' + str(i).zfill(9),
                                pos_list1, pos_list2, pos_list3, pos_list4, read1, read2)
