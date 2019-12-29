## 一、 实验目的

1. 熟悉并掌握fna、fastq、sam文件格式
2. 掌握自索引构建算法
3. 运用索引进行read精准匹配即编辑距离k以内的模糊匹配

## 二、   实验原理

### 0. fna和fastq文件读入

fna文件除了第一行外都是碱基序列，我们只需要把除了第一行的其他行所有字符连接到一起就可以得到参考序列。

fastq文件中包含多个碱基序列。fastq1和fastq2两个文件成对出现。对于同一个名字测试fastq，有一对pair-end测试序列。起哄每个测试序列只有第二行代表碱基序列。因此每隔4行有一个有用的碱基序列。这样读入即可。

 ### 1. 构造BWT索引

BWT，数据转换算法，其实也是一种压缩算法，基本思想就是找到字符串的重复部分来进行压缩，还可以用来进行序列比对。BWT会将字符串转换成一个类似的字符串，这个字符串相邻的字符相同的可能性很大，这样，我们就可以对数据进行压缩了。

在该实验中，我们构造BWT索引来进行序列比对而非压缩。

BWT编码压缩步骤如下：

1. 首先对要转换的字符串，添加一个不在字符串里的ASCII码表里最小的字符。如 AGGAGC ——> $AGGAGC，添加了 \$。
2. 对字符串进行依次循环移位，得到一系列的字符串，如果字符串长度为 n， 就可以得到n个字符串，如下面图里的第二列所示。
3. 对2中的位移后的一系列的字符串按照ASCII进行排序，如下图的第四列所示，第三列是排序后的字符串的原index位置。
4. 取位移后的一系列字符串的首字母出来作为 F 列， 最后一个字母作为 L 列。如下图 F 列 和L 列所示。
5. **L 列就是最后的编码结果。**

![1577451799109](C:\Users\Stuyxr\AppData\Roaming\Typora\typora-user-images\1577451799109.png)

由编码过程，参考上图，其实可以发现BWT编码有三个特性(循环位移决定)，

1. L 列的第一个元素是源字符串的最后一个元素。
2. 循环位移可知，同一行的 F 列和 L 列的元素在源字符串里是相邻的，而且 L 列元素的下一个字符就是 同行里 F 列的元素。
3. 同一种字符在 F 列和 L 列里的rank是一样的，比方说， F 列里的第二个 A 和 L 列里的第二个 A 在源字符串里是同一个A。 F 列里的第一个 G 和 L 列里的第一个 G 在源字符串里是同一个G, rank如下图所示。

![rank](https://chengjunwen.github.io/images/bwt/FLrank.png)

具体BWT索引的构造方法，我们借助后缀数组来实现。不难发现，设后缀数组为`sa[]`，原字符串为是`s`，那么`bwt[i] = s[sa[i]-1]`。那么我们只要构造出后缀数组就可以求出BWT索引了。

在该实验中，我后缀数组用的是SA-IS算法。SA-IS可以在线性时间$O(n)$内计算出后缀数组。又由于python跑该大肠杆菌参考基因组的速度实在捉急，我用c++实现的后缀数组部分，见附件`sais.cpp`。该程序将后缀数组输出到`data/sa.txt`和`data/rev_sa.txt`（翻转）中。接着python项目从`data/sa.txt`和`data/rev_sa.txt`中读取后缀数组，而不是用python计算，从而提高运行速度。

SA-IS算法的伪代码如下：

```
function SA-IS(S):
    t = bool[]
    S1 = int[]
    P = int[]
    bucket = int[]
    扫描倒序字符串确定每一个后缀的类型 -> t
    扫描t数组确定所有的 LMS 子串 -> P
    对所有的 LMS 子串进行诱导排序
    对每一个 LMS 子串重新命名，生成新的串 S1

    if S1 中的每一个字符都不一样:
        直接计算 SA1
    else
        SA1 = SA-IS(S1)  # 递归计算 SA1

    利用 SA1 来进行诱导排序，计算 SA
    return SA
```

SA-IS 算法是基于诱导排序这种思想。基本想法就是将问题的规模缩小，通过解决更小的问题，获取足够信息，就可以快速的解决原始问题。从这里也可以看出，这一过程需要递归处理子问题。

### 2. BWT索引解码出原字符串

有了BWT索引，我们考虑如何解码还原出原字符串，方法如下：

1. 由 L 列 得到 F 列。因为L 列 和F 列其实都是源串的字符的不同排列方式，但是我们知道 F 列是按照 ASCII码排序的，所以从 L 就可以推出 F 。
2. 根据第一个性质，我们可以得到源串的最后一个字符是 L 列的第一个字符，作为当前字符(下面依次往前递推)。
3. 依据上一步得到的作为当前字符， 根据第三个性质，我们可以得到同一个字符在 F 列中的位置，作为当前字符。
4. 依据 F 列里的当前字符，根据第二个性质，我们可以得到当前字符的上一个字符是同行里的 L 列里的元素,将新增字符作为当前字符，然后跳转到第 3 步。
5. 直到所有字符全部推算出来。

过程如下：

![jiema](https://chengjunwen.github.io/images/bwt/LFmapping.png)

解密部分伪代码如下：

![1577453076952](C:\Users\Stuyxr\AppData\Roaming\Typora\typora-user-images\1577453076952.png)

其中，C[c]的定义是，字典序小于字符c的所有字符个数。Occ[c,r]表示在BWT中第r行之前出现字符c的个数。

例如，下面这个例子，原串为acaacg\$，BWT为gc\$aaac，C[a]=0,C[c]=3,C[g]=5。Occ[a,5] =2。从0开始索引。

![1577453176547](C:\Users\Stuyxr\AppData\Roaming\Typora\typora-user-images\1577453176547.png)

按照上述代码执行还原过程：

![1577453334102](C:\Users\Stuyxr\AppData\Roaming\Typora\typora-user-images\1577453334102.png)

### 3. 字符串精准匹配FM-Index

算法流程如下：

![1577453377781](C:\Users\Stuyxr\AppData\Roaming\Typora\typora-user-images\1577453377781.png)

LFC（r,c）是对LF（r）的改进，在LF（r）中，第r行对应的只能是字符c，而LFC（r,c）中可以是任意的c。

我们将read从后向前匹配，在准确匹配过程中，我们维护一个区间[sp, ep)，该区间内的为当前匹配的位置的区间。下图是在acaacg中查找aac的实例。第一步我们匹配c，sp和ep分别是4和6；接着匹配ac，sp和ep由LFC函数计算而得，结果是sp=2，ep=4；接着匹配aac，sp=1，ep=2。

这里用的一个性质是如果S的匹配区间是[sp, ep)，那么xS的匹配范围也一定是一个区间[sp', ep')。

![1577453472707](C:\Users\Stuyxr\AppData\Roaming\Typora\typora-user-images\1577453472707.png)

### 4. BWA算法：编辑距离x以内的字符串匹配

在一个节点，测序序列的碱基与参考序列的碱基一定是几种情况之一：

1. 正确，与原序列一致
2. 错误，由另外三种碱基替换而成
3. 错误，该碱基前置一个或多个碱基被删除
4. 错误，该碱基是被插入的

四种情况下，第一种情况没有罚分，而剩下三种情况都发生了一次错误，都罚一分，并且需要参考序列或者测序序列向前移动一下：

![1577455883512](C:\Users\Stuyxr\AppData\Roaming\Typora\typora-user-images\1577455883512.png)

最朴素的做法就是搜索，枚举每一种错误，但是这样的话算法时间复杂度过高，无法接受。因此我们考虑剪枝：匹配的时候，每个位置设置一个错误上限$D_i$，当你的错误个数超过这个上限$D_i$时，我们就不接着向下搜索了，这样就可以大大降低时间复杂度。

那么我们如何计算$D_i$。我们计算$D_i$为参考序列和测试序列`P[:i]`之间不一样的字符的个数，那么这个$D_i$就是错误上限。$D_i$可以类似准确匹配的方法取匹配，知道匹配不到，从头再来，将错误上限+1。

计算错误上限$D$的伪代码如下：

![1577458546086](C:\Users\Stuyxr\AppData\Roaming\Typora\typora-user-images\1577458546086.png)

接着，进行搜索：InexRecur 函数的五个参数分别是：测序序列 W ；当前匹配位 i 最多允许错误数 z ；以及 BWT 位置上下限 k 、 l。

搜索部分首先对于$z<D_i$的情况剪枝，对于搜到头的部分$i<0$返回结果。接着枚举错误类型：插入错误、删除错误、碱基替换，以及没有错误的情况。）除了正确匹配不处罚之外，其他情况 $z-1$。）删除错误 i 不动，其他情况匹配位左移一位。伪代码如下图：

![1577458735377](C:\Users\Stuyxr\AppData\Roaming\Typora\typora-user-images\1577458735377.png)

BWA算法的主程序的伪代码如下：

![1577458865971](C:\Users\Stuyxr\AppData\Roaming\Typora\typora-user-images\1577458865971.png)

BWA算法思想是：

1. 对于查找序列上的每个位置，都使用 [A C G T] 以及插入删除共计 6 种情况进行查找。
2. 6 种情况中，对于序列一致的分支，不积累错误，其他 5 个分支各积累一个错误。
3. 错误达到上限，剪枝。

### 5. SAM文件输出

1. 第一列：read name，read的名字通常包括测序平台等信息；

2. 第二列：sum of flags，比对flag数字之和，比对flag用数字表示，分别为：

3. 第三列：RNAM，reference sequence name，实际上就是比对到参考序列上的染色体号。若是无法比对，则是*；

   - 0 read是Single-read且正向比对

   - 1 read是pair中的一条（read表示本条read，mate表示pair中的另一条read）
   - 2 pair一正一负完美的比对上
   - 4 这条read没有比对上
   - 8 mate没有比对上
   - 16 这条read反向比对
   - 32 mate反向比对
   - 64 这条read是read1
   - 128 这条read是read2
   - 256 第二次比对
   - 512 比对质量不合格
   - 1024 read是PCR或光学副本产生
   - 2048 辅助比对结果

   通过这个和可以直接推断出匹配的情况。假如说标记不是以上列举出的数字，比如说83=（64+16+2+1），就是这几种情况值和。

4. 第四列：position，read比对到参考序列上，第一个碱基所在的位置。若是无法比对，则是0；

5. 第五列：Mapping quality，比对的质量分数，越高说明该read比对到参考基因组上的位置越唯一；

6. 第六列：CIGAR值，read比对的具体情况：

   - “M”表示 match或 mismatch；
   - “I”表示 insert；
   - “D”表示 deletion；
   - “N”表示 skipped（跳过这段区域）；
   - “S”表示 soft clipping（被剪切的序列存在于序列中）；
   - “H”表示 hard clipping（被剪切的序列不存在于序列中）；
   - “P”表示 padding；
   - “=”表示 match；
   - “X”表示 mismatch（错配，位置是一一对应的）；

7. 第七列：MRNM(chr)，mate的reference sequence name，实际上就是mate比对到的染色体号，若是没有mate，则是*；

8. 第八列：mate position，mate比对到参考序列上的第一个碱基位置，若无mate,则为0；

9. 第九列：ISIZE，Inferred fragment size.详见Illumina中paired end sequencing 和 mate pair sequencing，是负数，推测应该是两条read之间的间隔(待查证)，若无mate则为0；

10. 第十列：Sequence，就是read的碱基序列，如果是比对到互补链上则是reverse completed eg. CGTTTCTGTGGGTGATGGGCCTGAGGGGCGTTCTCN

11. 第十一列：ASCII，read质量的ASCII编码。

12. 第十二列之后：Optional fields，可选的自定义区域

将对应的参数传入，输出文件

![1577459466771](C:\Users\Stuyxr\AppData\Roaming\Typora\typora-user-images\1577459466771.png)

## 三、测试结果及分析

 准确匹配：

![1577459495035](C:\Users\Stuyxr\AppData\Roaming\Typora\typora-user-images\1577459495035.png)

编辑距离12以内的匹配：

![1577459831416](C:\Users\Stuyxr\AppData\Roaming\Typora\typora-user-images\1577459831416.png)

分析：准确匹配时，有的测序片段能匹配上，有的测序片段不能匹配上，可能放生了基因突变；而容错下几乎所有序列都能匹配上，说明基因突变概率不大。

## 四、经验体会

 在本次实验中，我学会了对于在生物信息中对于基因序列的比对方法，实现了BWT索引、测序片段精准匹配和模糊匹配的方法。对于科技前沿技术有了进一步的认知。

## 五、  附录：源代码（带注释）

```c++
// Author: Stuyxr
// SA-IS算法构造后缀数组
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
```

```python
# Author: Stuyxr
# 加载和输出文件
def load_file(filename):
    '''
    从文件中读取数据
    :param filename: 文件名字
    :return: 读取出来的原始数据
    '''
    with open(filename, 'r') as fp:
        data = fp.read()
    return data


def load_ref(file):
    '''
    对参考序列进行一些处理，去掉描述信息并将数据全部转换为大写
    :param filename: 存储参考序列的文件
    :return: 全部大写的参考序列
    '''
    data = load_file(file)
    return data[data.index('\n') + 1:].replace('\n', '').upper()
    # [:100000]


def load_reads(file1, file2, num=1000):
    '''
    对需要比对的read进行处理，包括将文件中的read抽取出来，并将数据全部转换为大写。
    这里因为是pair-end的read，所以同时将两个文件中的reads分别取出来
    :param num: 读取read的个数
    :param file1: 存储read的文件
    :param file2: 存储read的文件
    :return: 提取出来的read列表
    '''
    with open(file1, 'r') as fp:
        data1 = fp.read()
    with open(file2, 'r') as fp:
        data2 = fp.read()
    data1 = data1.split('\n')
    data2 = data2.split('\n')
    reads1 = []
    reads2 = []
    for i in range(int(num)):
        reads1.append(data1[i * 4 + 1].upper())
        reads2.append(data2[i * 4 + 1].upper())
    return reads1, reads2


def calculate_flag(pos_list1, pos_list2, pos_list3, pos_list4, number):
    """
    计算flag值
    :param pos_list1: 测序片段1匹配到的位置list
    :param pos_list2: 测序片段2匹配到的位置list
    :param pos_list3: 测序片段3(1的互补反序序列)匹配到的位置list
    :param pos_list4: 测序片段4(2的互补反序序列)匹配到的位置list
    :param number: number=1时表示当前测试序列是测试序列1，否则当前测试序列是测试序列2
    :return: flag值
    """
    flag = 1
    if len(pos_list1) == 0:
        flag += 4
    if len(pos_list2) == 0:
        flag += 8
    if len(pos_list3) != 0:
        flag += 16
    if len(pos_list4) != 0:
        flag += 32
    if number == 1:
        flag += 64
    else:
        flag += 128
    return flag


def getline(qname, flag, pos, pnext, seq, mapq='60', cigar='=', rname='NC_008253', tlen='0', qual='*'):
    """
    格式化一行信息
    :param qname: 测序片段的名字
    :param flag: 标志位
    :param pos: 偏移位置
    :param pnext: 另一个pair-end序列的偏移位置
    :param seq: 基因序列
    :param mapq: 质量分数
    :param cigar: cigar值
    :param rname: 参考序列名字
    :param tlen: ISIZE值
    :param qual: read质量的ASCII编码
    :return: 格式化后的一行字符串
    """
    return qname + '\t' + str(flag) + '\t' + rname + '\t' + str(pos) + '\t' + mapq + '\t' + cigar + '\t' \
           + rname + '\t' + str(pnext) + '\t' + tlen + '\t' + seq + '\t' + qual + '\n'


def write_line(pos_list1, pos_list2, fp, read_name, flag, seq):
    """
    写入文件的一行
    :param pos_list1: 测序片段1匹配到的位置list
    :param pos_list2: 测序片段2匹配到的位置list
    :param fp: 文件
    :param read_name: 测序片段名字
    :param flag: 标志位
    :param seq: 测序片段
    :return: None
    """
    pnext = 0 if len(pos_list2) == 0 else pos_list2[0]
    if len(pos_list1) != 0:
        for pos in pos_list1:
            fp.write(getline(qname=read_name, flag=flag, pos=pos, pnext=pnext, seq=seq))
    else:
        fp.write(getline(qname=read_name, flag=flag, pos=0, pnext=pnext, seq=seq))


def write_head(file_name):
    """
    写入文件头部
    :param file_name: 文件名
    :return: None
    """
    with open(file_name, 'w') as fp:
        fp.write('@HD VN:1.6 SO:coordinate' + '\n' + '@SQ SN:ref LN:45' + '\n')


def write_sam(file_name, read_name, pos_list1, pos_list2,  pos_list3, pos_list4, read1, read2):
    """
    写入sam文件
    :param file_name: 文件名
    :param read_name: 测序片段名字
    :param pos_list1: 测序片段1匹配到的位置list
    :param pos_list2: 测序片段2匹配到的位置list
    :param pos_list3: 测序片段3(1的互补反序序列)匹配到的位置list
    :param pos_list4: 测序片段4(2的互补反序序列)匹配到的位置list
    :param read1: 测序片段1
    :param read2: 测序片段2
    :return: None
    """
    with open(file_name, 'a') as fp:
        flag1 = calculate_flag(pos_list1, pos_list2,  pos_list3, pos_list4, 1)
        flag2 = calculate_flag(pos_list2, pos_list1,  pos_list4, pos_list3, 2)
        write_line(pos_list1, pos_list2, fp, read_name, flag1, read1)
        write_line(pos_list2, pos_list1, fp, read_name, flag2, read2)
```

```python
# Author: Stuyxr
# BWA算法
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
```

