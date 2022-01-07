这个文件夹存放了目前实验跑的一些数据集，其中一些可能存在问题:
1. amazon-2008 directed
   @网页记录的数据  nodes=735323 arcs=5158388
   @实验跑的数据    nodes=735323 arcs=7046944 self\_loop=0
2. **cnr-2000** directed
   @网页记录的数据  nodes=325557 arcs=3216152
   @实验跑的数据    nodes=325557 arcs=5565380 **self\_loop=87442**
3. dblp-2010 undirected
    @网页记录的数据  nodes=326186 arcs=1615400
    @实验跑的数据    nodes=326186 arcs=1615400 self\_loop=0
4. dblp-2011 undirected
    @网页记录的数据  nodes=986324 arcs=6707236
    @实验跑的数据    nodes=986324 arcs=6707236 self\_loop=0
5. **enron** directed
    @网页记录的数据  nodes=69244 arcs=276143
    @实验跑的数据    nodes=69244 arcs=510433 **self\_loop=1535**
6. wordassociation-2011 directed
    @网页记录的数据  nodes=10617 arcs=5158388
    @实验跑的数据    nodes=10617 arcs=127576 self\_loop=0
7. **in-2004** directed
    @网页记录的数据  nodes=1382908 arcs=16917053
    @实验跑的数据    nodes=1382908 arcs=27560356 self\_loop=377410

加粗显示的数据集有**cnr-2000**，**enron**和**in-2004**三个，这三个都是有向图，但是在转换成无向图的过程中存在问题，即存在self-loop(自环边)，这应该是错误的。
