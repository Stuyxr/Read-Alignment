import sais
import loading

dict1 = {'$': 0, 'A': 1, 'C': 2, 'G': 3, 'T': 4, 'Z': 5}
dict2 = ['$', 'A', 'C', 'G', 'T', 'Z']
k = ord('T') + 1


def get_bwt_array(s):
    """
    C[c]的定义是，字典序小于字符c的所有字符个数。
    Occ[r, c]表示在bwt[]中第r行之前出现字符c的个数。
    :param s:
    :return:
    """
    bwt = []
    C = [0] * k
    suffix_arr = sais.sa_is(s, k)
    for i in range(len(s)):
        for x in dict1.keys():
            if s[i] < x:
                C[dict1[x]] += 1
    n = len(s)
    for i in range(n):
        bwt.append(s[(suffix_arr[i] - 1 + n) % n])
    Occ = []
    for i in range(len(bwt)):
        if i == 0:
            Occ.append([0, 0, 0, 0, 0])
        else:
            Occ.append(list(Occ[i - 1]))
        Occ[i][dict1[bwt[i]]] += 1
    Occ.append(list(Occ[len(bwt) - 1]))
    for i in range(len(Occ) - 1):
        Occ[i][dict1[bwt[i]]] -= 1
    return bwt, C, Occ, suffix_arr


def lfc(r, c):
    c = dict1[c]
    return C[c] + Occ[r][c]


def r_lfc(r, c):
    c = dict1[c]
    return r_C[c] + r_Occ[r][c]


def exact_match(P):
    """
    字符串精确匹配
    :param P: 字符串P
    :return:
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
    print(sp, ep)
    return sp < ep


def calculate_d(P):
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
    print(i, z, sp, ep, P[i])
    if i < 0:
        return set([x for x in range(sp, ep)])
    if z < d[i]:
        return set()
    s = set()
    s = s | inex_recur(P, i-1, z-1, sp, ep)
    spp, epp = sp, ep
    for c in dict2[1:5]:
        sp = lfc(sp, c)
        ep = lfc(ep, c)
        if sp < ep:
            s = s | inex_recur(P, i, z-1, sp, ep)
            if c == P[i]:
                s = s | inex_recur(P, i-1, z, sp, ep)
            elif P[i]:
                s = s | inex_recur(P, i-1, z-1, sp, ep)
        sp = spp
        ep = epp
    return s


def inexact_search(P, z):
    global d
    d = calculate_d(P)
    return inex_recur(P, len(P)-1, z, 0, len(ref))


if __name__ == '__main__':
    # ref = loading.load_ref('NC_008253.fna') + '$'
    # reads1, reads2 = loading.load_reads('NC_008253_1.fastq', 'NC_008253_2.fastq')
    ref = 'CCGCA' + '$'
    bwt, C, Occ, sa = get_bwt_array(ref)
    r_ref = ref[-2::-1] + '$'
    _, r_C, r_Occ, _ = get_bwt_array(r_ref)
    pos = inexact_search('CGC', 1)
    print(pos)
    print([sa[x] for x in pos])
