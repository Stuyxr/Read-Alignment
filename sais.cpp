#include <bits/stdc++.h>
using namespace std;
template<size_t size>
struct SuffixArray
{
    bool type[size*2];
    int bucket[size], bucket1[size];
    int sa[size], rk[size], ht[size];
    bool isLMS(const int i, const bool *type) { return i > 0 && type[i] && !type[i - 1]; }
    template<class T>
    void inducedSort(T s, int *sa, const int len, const int sigma, const int bucketSize, bool *type, int *bucket, int *cntbuf, int *p) {
		/** 
		 * 诱导排序(从*型诱导到L型、从L型诱导到S型)
		 * 调用之前应将*型按要求放入SA中
		 **/
        memset(bucket, 0, sizeof(int) * sigma);
        memset(sa, -1, sizeof(int) * len);
        for (int i = 0; i < len; i++) bucket[s[i]]++;
        cntbuf[0] = bucket[0];
        for (int i = 1; i < sigma; i++) cntbuf[i] = cntbuf[i - 1] + bucket[i];
        for (int i = bucketSize - 1; i >= 0; i--) sa[--cntbuf[s[p[i]]]] = p[i];
        for (int i = 1; i < sigma; i++) cntbuf[i] = cntbuf[i - 1] + bucket[i - 1];
        for (int i = 0; i < len; i++) if (sa[i] > 0 && !type[sa[i] - 1]) sa[cntbuf[s[sa[i] - 1]]++] = sa[i] - 1;
        cntbuf[0] = bucket[0];
        for (int i = 1; i < sigma; i++) cntbuf[i] = cntbuf[i - 1] + bucket[i];
        for (int i = len - 1; i >= 0; i--) if (sa[i] > 0 && type[sa[i] - 1]) sa[--cntbuf[s[sa[i] - 1]]] = sa[i] - 1;
    }
    template<typename T>
    inline void sais(T s, int *sa, int len, bool *type, int *bucket, int *bucket1, int sigma) {
		/** 
		 * SA-IS主体
		 * S是输入字符串，length是字符串的长度, SIGMA是字符集的大小
		 **/
        int  bucketSize = 0, cnt = 0, p = -1, x, *cntbuf = bucket + sigma;
        type[len - 1] = 1;
        for (int i = len - 2; i >= 0; i--) type[i] = s[i] < s[i + 1] || (s[i] == s[i + 1] && type[i + 1]);
        for (int i = 1; i < len; i++) if (type[i] && !type[i - 1]) bucket1[bucketSize++] = i;
        inducedSort(s, sa, len, sigma, bucketSize, type, bucket, cntbuf, bucket1);
        for (int i = bucketSize = 0; i < len; i++) if (isLMS(sa[i], type)) sa[bucketSize++] = sa[i];
        for (int i = bucketSize; i < len; i++) sa[i] = -1;
        for (int i = 0; i < bucketSize; i++) {
            x = sa[i];
            for (int j = 0; j < len; j++)
            {
                if (p == -1 || s[x + j] != s[p + j] || type[x + j] != type[p + j]) { cnt++, p = x; break; }
                else if (j > 0 && (isLMS(x + j, type) || isLMS(p + j, type))) break;
            }
            x = (~x & 1 ? x >> 1 : x - 1 >> 1), sa[bucketSize + x] = cnt - 1;
        }
        for (int i=len-1,j=len-1; i >= bucketSize; i--) if (sa[i] >= 0) sa[j--] = sa[i];
        int *s1 = sa + len - bucketSize, *bucket2 = bucket1 + bucketSize;
        if (cnt < bucketSize) sais(s1, sa, bucketSize, type + len, bucket, bucket1 + bucketSize, cnt);
        else for (int i = 0; i < bucketSize; i++) sa[s1[i]] = i;
        for (int i = 0; i < bucketSize; i++) bucket2[i] = bucket1[sa[i]];
        inducedSort(s, sa, len, sigma, bucketSize, type, bucket, cntbuf, bucket2);
    }
    inline void getHeight(const char *s, const int len, const int *sa) {
        for (int i = 0, k = 0; i < len; i++) {
            if (rk[i] == 0) k = 0;
            else {
                if (k > 0) k--;
                int j = sa[rk[i] - 1];
                while (i + k < len && j + k < len && s[i + k] == s[j + k]) k++;
            }
            ht[rk[i]] = k;
        }
    }
    template<class T>
    inline void init(T s, const int len, const int sigma) {
        sais(s, sa, len, type, bucket, bucket1, sigma);
        for (int i = 1; i < len; i++) rk[sa[i]] = i;
        getHeight(s, len, sa);
    }
};
const int MAXN = 1e7 + 10;
char s[MAXN], str[MAXN];
int len;
SuffixArray<MAXN>sa;
string sss;

int main()
{
	/**
	 * 输出到data/sa.txt或者data/rev_sa.txt中
	 **/
	freopen("data/NC_008253.fna", "r", stdin);
	freopen("data/rev_sa.txt", "w", stdout);
	getline(cin, sss);
	while(scanf("%s", str) != EOF) {
		strcat(s, str);
	}
	int n = strlen(s);
	for(int i = 0; i < n - 1 - i; i++)
		swap(s[i], s[n - 1 - i]);
	strcat(s, "$");
    len=strlen(s);
    sa.init(s,len +1,256);
    for(int i=1;i<=len;i++)
    {
        printf("%d%c",sa.sa[i],i==len?'\n':' ');
    }
    return 0;
}