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
