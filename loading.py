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
