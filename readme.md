# 基因本体相似度计算项目说明

该项目用于计算基因本体之间的相似度.并通过将基因术语的相似度累积为基因的相似度,最终通过生物路径计算对数差异倍数LFC进行评分.

## 项目结构
项目文件夹格式如下:
```
.
├── data                    //数据文件夹
│   └── ...
├── debug                   //debug文件夹,生成的二进制文件
│   └── GeneSimilarity
├── dir                     //其他数据文件夹
│   └── ...
├── include                 //头文件
│   └── ...
├── Makefile                //Makefile文件
├── obj                     //编译过程中间文件文件件
│   └── ...
├── readme.md               //说明文件
├── result                  //结果文件夹
│   └── ...         
├── src                     //源代码文件夹
│   └── ...
└── 注释说明.md              //注释文件(*.gaf)说明.
```

## 数据文件夹`data`
```
data/
├── child.buf
├── descendant.buf
├── ec.tab
├── gene.gaf
├── netmat.txt
├── net.txt
├── obo.buf
├── onto.obo
├── pair.txt
└── path.buf
```
上面的文件中,有四个文件是基础文件,其他的文件都是通过这四个文件生成的缓存文件.  
这四个文件分别是:
1. `ec.tab`  
    生物路径文件夹,存储基因与EC编号的数据,用于最终的LFC求值,来评估算法.
2. `gene.gaf`  
    基因注释文件,存储着基因和本体之间的注释关系,详细信息查看`注释说明.md`
3. `net.txt`
    基因功能网络数据,co-function的结果,这个文件是已有的全部的基因和连接.
4. `onto.obo`
    本体文件,保存着基因本体的信息.

**上面的基因注释文件和基因网络数据文件全部都是酿酒酵母的数据**

### 网络数据文件派生文件
通过`net.txt`派生出了一个矩阵,`netmat.txt`,将数据转化为一个矩阵.
**该矩阵是黄金矩阵,只包含了网络数据中可靠性最高的边,只使用了5800+基因中的4172个基因.这4000+基因的名称位于`dir/map.txt`中,按照行号索引.**

### 本体数据文件派生文件
通过`obo.net`中派生的文件有三个
1. `child.buf` 
    本体的直接子节点集合
2. `descendant.buf`
    本体的所有子孙节点集合
3. `path.buf`
    两个节点之间所有路径的节点
4. `obo.buf`
    本体结构缓存

**上述的文件全部都是为了加速计算,来避免反复的计算.**

### `pair.txt`
这个文件没有用处,来自于PJJ14年论文的数据.


## 中间数据文件夹`dir`
本文件夹中的数据都不太重要.

`result.mat`文件是使用0.5的重启参数随机游走跑出来的矩阵.算是一个遗留数据.  
`path.pair`是计算`data/path.buf`时的中间文件,用于确定哪一些路径需要计算,并去重复.
`map.txt`是基因名称和矩阵索引的映射文件,在从矩阵中读取网络数据的时候才需要使用.

## `result`数据文件夹
这个文件夹中的文件全部都是计算生成的文件,全部都可以重新计算出来.可有可无.

## 头文件文件夹`include`
```
include/
├── anno.h
├── defs.h
├── matrix.h
├── sim.h
└── term.h
```
目前只有五个文件,将来会将`sim.h`按照**Class**拆分开来,并在头文件中加上详细的注释和说明.

另外四个文件都是相关的数据结构定义.

## 源代码文件夹`src`
```
src/
├── gene.cpp        //计算基因相似度
├── lfc.cpp         //计算最终的lfc得分
├── main.cpp        //组合方法进行计算
├── term.cpp        //计算术语相似度(通过hashmap获取网络数据)
└── term_index.cpp  //计算术语相似度(通过矩阵的索引来获取网络数据)
```
文件夹中的文件按照计算的步骤分开放置.
按照LFC,基因相似度,术语相似度,来放置.
**为了加速计算,去除掉了基因相似度中的除一法部分,理论速度能够提升33倍.术语相似度计算部分添加了hash改数组的搜索,速度提升也极为明显.**

**矩阵填充部分的代码暂时没有编写,下一步将会进行补充.**

详细的计算步骤参见论文.



