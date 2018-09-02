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
