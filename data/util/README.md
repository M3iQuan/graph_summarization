这个文件夹里存放了下载数据集和使用数据集的方法
1. 首先使用 download.sh 脚本下载需要的数据集，注意替换 basename 变量;
2. 下载完成后使用 check.sh 脚本检查数据是否正常，以免下载时出错或数据发生来变化;
3. 然后使用 generate.sh 脚本生成所需的文件 basename.obl 和 basename.offsets，这两个文件用于在使用过程中能够实现random access，另外还生成了对称文件 \*.sym，主要是转换有向图成无向图。
