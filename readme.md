## 9.15 更新
在通过术语相似度计算基因相似度的时候,暂时忽略掉除一法,能够提升40倍左右的速度,精度影响暂时不明.

## 9.13 更新
通过使用二维数组来替代hashmap的网络数据能够在目前的基础上提高7倍左右的速度。

准备更换架构，数据不再集中管理，采用分散管理。

### 文件说明
data/                   数据文件夹，存储各类原始数据  
├── ec.tab                  生物路径数据  
├── gene.gaf                基因的注释数据  
├── net.txt                 基因功能网络的数据  
├── onto.obo                基因本体数据  
└── pair.txt                300对基因数据，pjj计算出来的现有数据，没有太大的参考价值，使用的各个文件都比PJJ的要新。

buf                     缓存数据文件夹，用于存放各类缓存数据，用于加速计算  
├── child.txt               节点到所有子孙节点  
└── son.txt                 节点到所有直系子节点  



## 9.8更新
添加LOG类，实时跟踪程序运行的位置。


# 说明
基因本体相似度计算的CPP版本。
在Java版本上进行优化升级，将各个模块之间解耦。使得各个模块都是独立的部分，以便于后期的改进。
主体上分成两部分，一部分是数据存储以及读取。另一部分负责计算。

## 数据存取
1. 读取文件
    读取对应的obo文件，gaf文件，net文件，ec文件，并转换为对应的数据结构
2. 缓存数据
    将中间数据尽可能的缓存，尤其是本体图的结构。等以加速计算。
3. 共享数据
    将数据使用const static对象保存起来。

**用于搜索的数据尽量使用无序hashmap。来加速搜索**


## 计算

1. 矩阵运算

2. 基因集合功能距离$D(G_a,G_b)$

3. 路径注释信息$f(t_a,t_b,p)$和$U(t_a,t_b,p)$



4. LCA特异性$h(t_a,t_b)$

5. 公共祖先集合`P`

6. 基因相似度$GS(g_1,g_2)$

7. 术语与术语集合集合相似度$Sim(t1,T_y)$
 $= max_{t_y \in  T_y}Sim_{term}(t,t_y)$

8. LFC累加
