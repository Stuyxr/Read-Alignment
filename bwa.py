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
    return bwt, C, Occ


def lfc(r, c):
    c = dict1[c]
    return C[c] + Occ[r][c]


def exact_match(P):
    """
    字符串精确匹配
    :param P: 字符串P
    :return:
    """
    p = len(P) - 1
    c = dict1[P[p]]
    sp = C[c]
    ep = C[c + 1]
    i = p - 1
    while sp < ep and i >= 0:
        c = P[i]
        sp = lfc(sp, c)
        ep = lfc(ep, c)
        i = i - 1
    print(sp, ep)
    return sp < ep


def calculate_d(W):
    print(r_ref)
    k = 0
    l = len(r_ref)
    z = 0
    d = [0] * len(W)
    for i in range(len(W)):
        c = dict1[W[i]]
        k = C[c] + rOcc[k][c]
        l = C[c] + rOcc[l][c]
        if k >= l:
            k = 0
            l = len(r_ref)
            z += 1
        d[i] = z
    return d


if __name__ == '__main__':
    # ref = loading.load_ref('NC_008253.fna') + '$'
    # reads1, reads2 = loading.load_reads('NC_008253_1.fastq', 'NC_008253_2.fastq')
    ref = 'AGACAGA' + '$'
    bwt, C, Occ = get_bwt_array(ref)
    # print(exact_match('CTTG'))

    r_ref = ref[-2::-1] + '$'
    r_bwt, rC, rOcc = get_bwt_array(r_ref)
    print(calculate_d('CGC'))
