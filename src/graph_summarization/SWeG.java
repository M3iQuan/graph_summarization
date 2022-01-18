package graph_summarization;

import java.util.*;
import java.util.concurrent.ThreadLocalRandom;

public class SWeG extends  Summary{

    // 用于对顶点的重新编号
    int[] h;
    // shingle数组，在顶点分组时使用
    int[] F;
    // 排序后的组别数组，满足 F[G[i]] <= F[G[i+1]] 即相同F值的在附近安排在同一个组内
    Integer[] G;
    // 找到G数组里的第一个组的第一个顶点的index
    int g_start;
    // 分组后的组别数量
    int num_groups;
    // 用于分组的二维数组，分别记录每个组的shingle值和第几组
    int[][] group_prop;

    // 用于计算合并的成功率
    ArrayList<Double> merged_success;
    HashMap<Integer, ArrayList<Integer>> total_groups_size;

    int signatureLength;

    /**
     * 构造函数，继承了Summary父类的构造函数，同时初始化自己的数据结构
     *
     * @param basename 数据集的基本名字
     * @throws Exception
     */
    public SWeG(String basename) throws Exception {
        super(basename);
        merged_success = new ArrayList<>();
        total_groups_size = new HashMap<>();
    }

    /**
     * 构造函数，继承了Summary父类的构造函数，同时初始化自己的数据结构
     *
     * @param basename 数据集的基本名字
     * @throws Exception
     */
    public SWeG(String basename, int signatureLength) throws Exception {
        super(basename);
        this.signatureLength = signatureLength;
        merged_success = new ArrayList<>();
        total_groups_size = new HashMap<>();
    }

    /**
     * 对顶点的一次重新编号 h: |V| -> |V|
     */
    private void randomPermutation() {
        h = new int[n];
        for (int i = 0; i < n; i++) {
            h[i] = i;
        }
        Random rnd = new Random();
        rnd = ThreadLocalRandom.current();
        for (int i = h.length - 1; i > 0; i--) {
            int index = rnd.nextInt(i + 1);
            int a = h[index];
            h[index] = h[i];
            h[i] = a;
        }
    }

    /**
     * 返回顶点u的shingle值, 方法是计算顶点集合 {u Union N(u)} 的最小 shingle 值
     *
     * @param u 顶点的编号
     * @return
     */
    private int shingleValue(int u) {
        int f_u = h[u];
//        int[] neighbors = Gr.successorArray(u);
        int[] neighbors = neighbors_[u];
        for (int v : neighbors) {
            if (f_u > h[v]) {
                f_u = h[v];
            }
        }
        return f_u;
    }

    /**
     * 传入参数k和随机生成的数组D, 为超点A计算其LSH签名, 使用的是单枚举算法
     * @param k LSH签名的长度
     * @param A 超点A的编号
     * @param rot_direction 随机生成的数组D
     * @return 返回超点A的LSH签名对象，为OnePermHashSig
     */
    public OnePermHashSig generateSignature(int k, int A, int[] rot_direction){
        // 通过参数k,确定签名的块数k_bins 和 每一块的长度 bin_size
        int k_bins = k;
        int bin_size = n / k_bins;
        if (n % k_bins != 0) { k_bins = k_bins + 1; }
        OnePermHashSig hashSig = new OnePermHashSig(k_bins);

        // A不是一个超点
        if (I[A] == -1) return hashSig;
        // 遍历超点的每个顶点v
        for (int v = I[A]; ; v=J[v]) {
            int[] neighbours = neighbors_[v];
//            int[] neighbours = Gr.successorArray(v);
            // 遍历顶点v的每个邻居j
            for (int neighbour : neighbours) {
                // 得到邻居j重排后的id
                int permuted_h = h[neighbour];
                // 以及在第几个bin
                int permuted_bin = permuted_h / bin_size;
                // 如果b_i==-1(即还没设置)或者比b_i还小的索引 (permuted_h%bin_size) 存在顶点，则重新设置b_i
                if (hashSig.sig[permuted_bin] == -1 || permuted_h % bin_size < hashSig.sig[permuted_bin]) {
                    hashSig.sig[permuted_bin] = permuted_h % bin_size;
                }
            }
            // 遍历完超点A的所有顶点就直接退出
            if(J[v]==-1) break;
        }

        // 开始处理块b_i是empty的情况，即rotation
        for (int A_bin = 0; A_bin < k_bins; A_bin++) {
            int direction = rot_direction[A_bin];
            // 如果b_i还没设置, 则利用rot_direction来确定往哪边采样
            if (hashSig.sig[A_bin] == -1) {
                int i = (A_bin + direction) % k_bins;
                if (i < 0) { i += k_bins; }
                int counter = 0;
                while (hashSig.sig[i] == -1 && counter < k_bins) {
                    i = (i + direction) % k_bins;
                    if (i < 0) { i += k_bins; }
                    counter++;
                }
                hashSig.sig[A_bin] = hashSig.sig[i];
            }
        }

        return hashSig;
    }

    /**
     * 计算一个组包含多少的顶点，即组的大小
     *
     * @param arr    排序后各组的shingle值，按照shingle值升序排序
     * @param key    要计算大小的那组的shingle值
     * @param st_pos 要计算大小的那组的开始index
     * @return
     */
    private int groupLength(int[] arr, int key, int st_pos) {
        int counter = 1;
        while (arr[st_pos++] == key && st_pos < arr.length - 1) counter++;
        return counter;
    }

    /**
     * 分组阶段，SWeG通过 shingle值对顶点进行划分
     * 分组完成后，可以通过遍历 G 得到每个组，相同组的G[i]值相等
     */
    @Override
    public double dividePhase() {
        logger_.info("开始分组, 采用的是Shingle方法");
        long startTime = System.currentTimeMillis();
        // 首先对顶点进行一个编号重排
        randomPermutation();
        // 初始化F数组, 用于存储每个顶点的shingle值
        F = new int[n];
        for (int A = 0; A < n; A++)
            F[A] = -1;
        for (int A = 0; A < n; A++) {
            // A不是一个超点
            if (I[A] == -1)
                continue;
            // 将超点A的shingle值先初始化成最大值，然后再逐渐通过超点包含的所有顶点的shingle值逐渐下降
            F[A] = n;
            for (int v = I[A]; ; v = J[v]) {
                int fv = shingleValue(v);
                if (F[A] > fv)
                    F[A] = fv;
                if (J[v] == -1)
                    break;
            }
        }

        // 对分组进行排序
        G = new Integer[n];
        for (int i = 0; i < n; i++) G[i] = i;
        Arrays.sort(G, (o1, o2) -> Integer.compare(F[o1], F[o2]));

        // 找到第一个组的开始index
        g_start = 0;
        while (F[G[g_start]] == -1)
            g_start++;

        // 计算组的数量
        num_groups = 0;
        int g = -1;
        for (int i = g_start; i < n; i++) {
            if (F[G[i]] != g) {
                num_groups++;
                g = F[G[i]];
            }
        }

        group_prop = new int[num_groups][2];
        g = -1;
        int counter = 0;
        for (int i = g_start; i < n; i++) {
            if (F[G[i]] != g) {
                group_prop[counter][0] = F[G[i]];
                group_prop[counter][1] = i;
                counter++;
                g = F[G[i]];
            }
        }

        return (System.currentTimeMillis() - startTime) / 1000.0;
    }

    /**
     * 合并阶段，SWeG算法在每个小组内的顶点采用Random方式进行合并
     */
    @Override
    public double mergePhase(double threshold) {
        System.out.println("# Merge Phase");
        System.out.println(String.format("Threshold=%5f", threshold));
        long startTime = System.currentTimeMillis();
        int idx = 0;
        // temp <- F[G]
        int[] temp = new int[n];
        for (int i = 0; i < n; i++) temp[i] = F[G[i]];

        // 开始遍历每个组
        for (int i = 0; i < num_groups; i++) {
            int st_position = group_prop[i][1];
            int group_size = groupLength(temp, group_prop[i][0], st_position) - 1;
            // 如果一个组只有一个顶点则直接跳过该组
            if (group_size < 2) continue;

            // Q是当前要进行合并的的组
            int[] Q = new int[group_size];
            int counter = 0;
            // 找到在Q组里面的所有超点
            for (int j = st_position; j < (st_position + group_size); j++) {
                Q[counter++] = G[j];
            }

            HashMap<Integer, HashMap<Integer, Integer>> hm = createW(Q, group_size);
            int initial_size = hm.size();
            while (hm.size() > 1) {
                Random rand = new Random();
                // 从组内随机找到一个超点A
                int A = rand.nextInt(initial_size);
                if (hm.get(A) == null)
                    continue;

                double max = 0;
                idx = -1;
                // 遍历组内其他顶点，找到与A的Jaccard Similarity最大的那个顶点
                for (int j = 0; j < initial_size; j++) {
                    if (hm.get(j) == null)
                        continue;
                    if (j == A) continue;
                    double jaccard_similarity = computeJacSim(hm.get(A), hm.get(j));
                    if (jaccard_similarity > max) {
                        max = jaccard_similarity;
                        idx = j;
                    }
                }
                if (idx == -1) {
                    hm.remove(A);
                    continue;
                }

                // 这里做了一个交换，目的是把编号较大的顶点合并到编号较小的顶点里面
                if (Q[A] > Q[idx]) {
                    int t = A;
                    A = idx;
                    idx = t;
                }

                // 计算两个顶点之间的合并收益
                double savings = computeSaving(hm.get(A), hm.get(idx), Q[A], Q[idx]);
                if (savings >= threshold) {
                    HashMap<Integer, Integer> w_update = updateW(hm.get(A), hm.get(idx));
                    hm.replace(A, w_update);
                    hm.remove(idx);
                    updateSuperNode(Q[A], Q[idx]);
                } else {
                    hm.remove(A);
                }
            }
        }
        return (System.currentTimeMillis() - startTime) / 1000.0;
    }

    /**
     * 测试函数，用于验证一下想法：
     *     1.每个组内的合并顶点数量 / 每个小组的顶点数量
     */
    public double mergePhase_test(int iteration, double threshold) {
        logger_.info(String.format("开始合并, 当前顶点的合并阈值:%5f", threshold));
        long startTime = System.currentTimeMillis();
        int idx = 0;
        // temp <- F[G]
        int[] temp = new int[n];
        for (int i = 0; i < n; i++) temp[i] = F[G[i]];

        HashMap<Integer, ArrayList<Integer>> merged_ = new HashMap<>();

        // 记录迭代次数i的所有size>2的组数量
        if(!total_groups_size.containsKey(iteration))
            total_groups_size.put(iteration, new ArrayList<>());

        // 开始遍历每个组
        for (int i = 0; i < num_groups; i++) {
            int st_position = group_prop[i][1];
            int group_size = groupLength(temp, group_prop[i][0], st_position) - 1;
            // 如果一个组只有一个顶点则直接跳过该组
            if (group_size < 2) continue;

            // group_size >= 2时，记录下group_size
            total_groups_size.get(iteration).add(group_size);

            // 能成功合并的对数
            int merged_pairs = 0;
            if (!merged_.containsKey(group_size)) {
                merged_.put(group_size, new ArrayList<>());
            }

            // Q是当前要进行合并的的组
            int[] Q = new int[group_size];
            int counter = 0;
            // 找到在Q组里面的所有超点
            for (int j = st_position; j < (st_position + group_size); j++) {
                Q[counter++] = G[j];
            }

            HashMap<Integer, HashMap<Integer, Integer>> hm = createW(Q, group_size);
            int initial_size = hm.size();
            while (hm.size() > 1) {
                Random rand = new Random();
                // 从组内随机找到一个超点A
                int A = rand.nextInt(initial_size);
                if (hm.get(A) == null)
                    continue;

                double max = 0;
                idx = -1;
                // 遍历组内其他顶点，找到与A的Jaccard Similarity最大的那个顶点
                for (int j = 0; j < initial_size; j++) {
                    if (hm.get(j) == null)
                        continue;
                    if (j == A) continue;
                    double jaccard_similarity = computeJacSim(hm.get(A), hm.get(j));
                    if (jaccard_similarity > max) {
                        max = jaccard_similarity;
                        idx = j;
                    }
                }
                if (idx == -1) {
                    hm.remove(A);
                    continue;
                }

                // 这里做了一个交换，目的是把编号较大的顶点合并到编号较小的顶点里面
                if (Q[A] > Q[idx]) {
                    int t = A;
                    A = idx;
                    idx = t;
                }

                // 计算两个顶点之间的合并收益
                double savings = computeSaving(hm.get(A), hm.get(idx), Q[A], Q[idx]);
                if (savings >= threshold) {
                    HashMap<Integer, Integer> w_update = updateW(hm.get(A), hm.get(idx));
                    hm.replace(A, w_update);
                    hm.remove(idx);
                    updateSuperNode(Q[A], Q[idx]);
                    merged_pairs++;
                } else {
                    hm.remove(A);
                }
            }
            merged_.get(group_size).add(merged_pairs);
        }

        double mm = 0.0;
        int nn = 0;
        for (Integer length : merged_.keySet()) {
            int sum = 0;
            int num = 0;
            for (Integer w : merged_.get(length)) {
                sum += w;
                num++;
            }
            mm += sum/(num*length*1.0);
            nn++;
        }
        merged_success.add(mm/nn);

        return (System.currentTimeMillis() - startTime) / 1000.0;
    }

    /**
     * 原始的分组方案, 不进行任何修改
     * @param iteration 迭代次数
     * @param print_iteration_offset 每隔多少次迭代就进行Encode一次查看压缩率
     */
    public void originTest(int iteration, int print_iteration_offset){
        long startTime = System.currentTimeMillis();
        Map<Integer, Map<Integer, Record>> it2Records = new HashMap<>();
        for (int it = 1; it <= iteration; it++) {
            List<List<Integer>> groups = new ArrayList<>();
            logger_.info(String.format("迭代轮数: %d", it));
            // DividePhase =============================================================================================
            {
                logger_.info("开始分组, 采用的是Local Sensitive Hash方法");
                long divideStartTime = System.currentTimeMillis();
                // 首先对顶点进行一个编号重排
                randomPermutation();
                // 初始化F数组, 用于存储每个顶点的shingle值
                int[] F_temp = new int[n];
                for (int A = 0; A < n; A++)
                    F_temp[A] = -1;
                for (int A = 0; A < n; A++) {
                    // A不是一个超点
                    if (I[A] == -1)
                        continue;
                    // 将超点A的shingle值先初始化成最大值，然后再逐渐通过超点包含的所有顶点的shingle值逐渐下降
                    F_temp[A] = n;
                    for (int v = I[A]; ; v = J[v]) {
                        int fv = shingleValue(v);
                        if (F_temp[A] > fv)
                            F_temp[A] = fv;
                        if (J[v] == -1)
                            break;
                    }
                }
                // 对分组进行排序
                Integer[] G_temp = new Integer[n];
                for (int i = 0; i < n; i++) G_temp[i] = i;
                Arrays.sort(G_temp, (o1, o2) -> Integer.compare(F_temp[o1], F_temp[o2]));
                // 找到第一个组的开始index
                int g_start_temp = 0;
                while (F_temp[G_temp[g_start_temp]] == -1)
                    g_start_temp++;

                int g = -1;
                ArrayList<Integer> Q = new ArrayList<>();
                for (int i = g_start_temp; i < n; i++) {
                    if (F_temp[G_temp[i]] != g) {
                       if(Q.size() > 1) groups.add(Q);
                       g = F_temp[G_temp[i]];
                       Q = new ArrayList<>();
                    }
                    Q.add(G_temp[i]);
                }
                if(Q.size() > 1) groups.add(Q);
                logger_.info(String.format("分组结束, 花费时间 %5f seconds", (System.currentTimeMillis() - divideStartTime) / 1000.0));
            }
            // End =====================================================================================================

            // MergePhase ==============================================================================================
            {
                double threshold = 1 / ((it + 1) * 1.0);
                logger_.info(String.format("开始合并, 当前顶点的合并阈值:%5f", threshold));
                // 记录mergePhase的运行时间
                long mergeStartTime = System.currentTimeMillis();
                // 记录迭代次数i的所有size>2的组数量
                if(!total_groups_size.containsKey(it))
                    total_groups_size.put(it, new ArrayList<>());
                Map<Integer, Record> recordMap = new HashMap<>();
                // 开始遍历每个组
                for (List<Integer> group : groups) {
                    long start = System.currentTimeMillis();
                    int success = 0;
                    // group_size >= 2时，记录下group_size
                    total_groups_size.get(it).add(group.size());
                    HashMap<Integer, HashMap<Integer, Integer>> hm = createW(group);
                    int initial_size = hm.size();
                    while (hm.size() > 1) {
                        Random rand = new Random();
                        // 从组内随机找到一个超点A
                        int A = rand.nextInt(initial_size);
                        if (hm.get(A) == null) continue;
                        // 变量max记录最大的Jaccard_similarity
                        double max = 0;
                        // 变量idx记录最大Jaccard_similarity的那个顶点
                        int idx = -1;
                        // 遍历组内其他顶点，找到与A的Jaccard Similarity最大的那个顶点
                        for (int j = 0; j < initial_size; j++) {
                            if (hm.get(j) == null)
                                continue;
                            if (j == A) continue;
                            double jaccard_similarity = computeJacSim(hm.get(A), hm.get(j));
                            if (jaccard_similarity > max) {
                                max = jaccard_similarity;
                                idx = j;
                            }
                        }
                        if (idx == -1) {
                            hm.remove(A);
                            continue;
                        }
                        // 这里做了一个交换，目的是把编号较大的顶点合并到编号较小的顶点里面
                        if (group.get(A) > group.get(idx)) {
                            int t = A;
                            A = idx;
                            idx = t;
                        }
                        // 计算两个顶点之间的合并收益
                        double savings = computeSaving(hm.get(A), hm.get(idx), group.get(A), group.get(idx));
                        if (savings >= threshold) {
                            HashMap<Integer, Integer> w_update = updateW(hm.get(A), hm.get(idx));
                            hm.replace(A, w_update);
                            hm.remove(idx);
                            updateSuperNode(group.get(A), group.get(idx));
                            success += 1;
                        } else {
                            hm.remove(A);
                        }
                    }
                    Double time = (System.currentTimeMillis() - start) / 1000.0;
                    if (!recordMap.containsKey(group.size())) {
                        Record r = new Record(group.size(), 1, success, time);
                        recordMap.put(group.size(), r);
                    } else {
                        Record r = recordMap.get(group.size());
                        r.add(1, success, time);
                    }

                }
                it2Records.put(it, recordMap);
                logger_.info(String.format("合并结束, 花费时间 %5f seconds", (System.currentTimeMillis() - mergeStartTime) / 1000.0));
            }
            // End =====================================================================================================
            if (it % print_iteration_offset == 0) {
                logger_.info(String.format("编码结束, 花费时间 %5f seconds", encodePhase_test()));
                evaluatePhase();
                logger_.info(String.format("到这里, 花费时间 %5f seconds", (System.currentTimeMillis() - startTime) / 1000.0));
            }
        }
        for (Integer it : total_groups_size.keySet()) {
            int max = 0;
            for (Integer length : total_groups_size.get(it)) {
                if (max <= length) {
                    max = length;
                }
            }
            logger_.info(String.format("迭代次数 %d: 小组最大顶点数量=%d", it, max));
        }
        logger_.info(String.format("程序运行结束, 总花费时间 %5f seconds", (System.currentTimeMillis() - startTime) / 1000.0));
    }

    /**
     * 对组内顶点数量大于maxGroupSize的小组进行顺序划分
     * @param group 组内数量大于maxGroupSize的小组
     * @param maxGroupSize 最大组内数量阈值
     * @return 按阈值切分后的小组
     */
    public List<List<Integer>> sequentialSplitGroup(List<Integer> group, int maxGroupSize){
        List<List<Integer>> groups = new ArrayList<>();
        int group_size = group.size();
        int num = group_size / maxGroupSize;
        if(group_size % maxGroupSize != 0) num += 1;
        for (int j = 0; j < num; j++) {
            int start = j * maxGroupSize;
            int end = Math.min((j + 1) * maxGroupSize, group_size);
            List<Integer> subGroup = new ArrayList<>();
            for (int k = start; k < end; k++) {
                subGroup.add(group.get(k));
            }
            groups.add(subGroup);
        }
        return  groups;
    }

    /**
     * 顺序切割的分组方案
     * @param iteration 迭代次数
     * @param print_iteration_offset 每隔多少次迭代就进行Encode一次查看压缩率
     */
    public void sequentialTest(int iteration, int print_iteration_offset, int maxGroupSize){
        long startTime = System.currentTimeMillis();
        Map<Integer, Map<Integer, Record>> it2Records = new HashMap<>();
        for (int it = 1; it <= iteration; it++) {
            List<List<Integer>> groups = new ArrayList<>();
            logger_.info(String.format("迭代轮数: %d", it));
            // DividePhase =============================================================================================
            {
                logger_.info("开始分组, 采用的是Local Sensitive Hash方法");
                long divideStartTime = System.currentTimeMillis();
                // 首先对顶点进行一个编号重排
                randomPermutation();
                // 初始化F数组, 用于存储每个顶点的shingle值
                int[] F_temp = new int[n];
                for (int A = 0; A < n; A++)
                    F_temp[A] = -1;
                for (int A = 0; A < n; A++) {
                    // A不是一个超点
                    if (I[A] == -1)
                        continue;
                    // 将超点A的shingle值先初始化成最大值，然后再逐渐通过超点包含的所有顶点的shingle值逐渐下降
                    F_temp[A] = n;
                    for (int v = I[A]; ; v = J[v]) {
                        int fv = shingleValue(v);
                        if (F_temp[A] > fv)
                            F_temp[A] = fv;
                        if (J[v] == -1)
                            break;
                    }
                }
                // 对分组进行排序
                Integer[] G_temp = new Integer[n];
                for (int i = 0; i < n; i++) G_temp[i] = i;
                Arrays.sort(G_temp, (o1, o2) -> Integer.compare(F_temp[o1], F_temp[o2]));
                // 找到第一个组的开始index
                int g_start_temp = 0;
                while (F_temp[G_temp[g_start_temp]] == -1)
                    g_start_temp++;

                int g = -1;
                ArrayList<Integer> Q = new ArrayList<>();
                for (int i = g_start_temp; i < n; i++) {
                    if (F_temp[G_temp[i]] != g) {
                        if(Q.size() > 1){
                            if(Q.size() <= maxGroupSize){
                                groups.add(Q);
                            } else {
                                logger_.info(String.format("当前小组顶点数量为%d, 大于阈值%d, 进行顺序切割", Q.size(), maxGroupSize));
                                List<List<Integer>> subGroups = sequentialSplitGroup(Q, maxGroupSize);
                                for (List<Integer> subGroup : subGroups) {
                                    if(subGroup.size() > 1) groups.add(subGroup);
                                }
                            }
                        }
                        g = F_temp[G_temp[i]];
                        Q = new ArrayList<>();
                    }
                    Q.add(G_temp[i]);
                }
                if(Q.size() > 1) groups.add(Q);
                logger_.info(String.format("分组结束, 花费时间 %5f seconds", (System.currentTimeMillis() - divideStartTime) / 1000.0));
            }
            // End =====================================================================================================

            // MergePhase ==============================================================================================
            {
                double threshold = 1 / ((it + 1) * 1.0);
                logger_.info(String.format("开始合并, 当前顶点的合并阈值:%5f", threshold));
                // 记录mergePhase的运行时间
                long mergeStartTime = System.currentTimeMillis();
                // 记录迭代次数i的所有size>2的组数量
                if(!total_groups_size.containsKey(it))
                    total_groups_size.put(it, new ArrayList<>());
                Map<Integer, Record> recordMap = new HashMap<>();
                // 开始遍历每个组
                for (List<Integer> group : groups) {
                    long start = System.currentTimeMillis();
                    int success = 0;
                    // group_size >= 2时，记录下group_size
                    total_groups_size.get(it).add(group.size());
                    HashMap<Integer, HashMap<Integer, Integer>> hm = createW(group);
                    int initial_size = hm.size();
                    while (hm.size() > 1) {
                        Random rand = new Random();
                        // 从组内随机找到一个超点A
                        int A = rand.nextInt(initial_size);
                        if (hm.get(A) == null) continue;
                        // 变量max记录最大的Jaccard_similarity
                        double max = 0;
                        // 变量idx记录最大Jaccard_similarity的那个顶点
                        int idx = -1;
                        // 遍历组内其他顶点，找到与A的Jaccard Similarity最大的那个顶点
                        for (int j = 0; j < initial_size; j++) {
                            if (hm.get(j) == null)
                                continue;
                            if (j == A) continue;
                            double jaccard_similarity = computeJacSim(hm.get(A), hm.get(j));
                            if (jaccard_similarity > max) {
                                max = jaccard_similarity;
                                idx = j;
                            }
                        }
                        if (idx == -1) {
                            hm.remove(A);
                            continue;
                        }
                        // 这里做了一个交换，目的是把编号较大的顶点合并到编号较小的顶点里面
                        if (group.get(A) > group.get(idx)) {
                            int t = A;
                            A = idx;
                            idx = t;
                        }
                        // 计算两个顶点之间的合并收益
                        double savings = computeSaving(hm.get(A), hm.get(idx), group.get(A), group.get(idx));
                        if (savings >= threshold) {
                            HashMap<Integer, Integer> w_update = updateW(hm.get(A), hm.get(idx));
                            hm.replace(A, w_update);
                            hm.remove(idx);
                            updateSuperNode(group.get(A), group.get(idx));
                            success += 1;
                        } else {
                            hm.remove(A);
                        }
                    }
                    Double time = (System.currentTimeMillis() - start) / 1000.0;
                    if (!recordMap.containsKey(group.size())) {
                        Record r = new Record(group.size(), 1, success, time);
                        recordMap.put(group.size(), r);
                    } else {
                        Record r = recordMap.get(group.size());
                        r.add(1, success, time);
                    }

                }
                it2Records.put(it, recordMap);
                logger_.info(String.format("合并结束, 花费时间 %5f seconds", (System.currentTimeMillis() - mergeStartTime) / 1000.0));
            }
            // End =====================================================================================================
            if (it % print_iteration_offset == 0) {
                logger_.info(String.format("编码结束, 花费时间 %5f seconds", encodePhase_test()));
                evaluatePhase();
                logger_.info(String.format("到这里, 花费时间 %5f seconds", (System.currentTimeMillis() - startTime) / 1000.0));
            }
        }
        for (Integer it : total_groups_size.keySet()) {
            int max = 0;
            for (Integer length : total_groups_size.get(it)) {
                if (max <= length) {
                    max = length;
                }
            }
            logger_.info(String.format("迭代次数 %d: 小组最大顶点数量=%d", it, max));
        }
        logger_.info(String.format("程序运行结束, 总花费时间 %5f seconds", (System.currentTimeMillis() - startTime) / 1000.0));
    }


    /**
     * 对于某个小组内的顶点进行分组
     * @param group 需要再进行划分的小组
     * @return 返回生成的组，每一个ArrayList都是一组，组内元素是超点的编号
     */
    public List<List<Integer>> hierarchicalSplitGroup(List<Integer> group, String method){
        int group_size = group.size();
        logger_.info(String.format("采用%s方式进行层次分组, 当前小组顶点数量为%d", method, group_size));

        List<List<Integer>> groups = new ArrayList<>();
        // 接着对顶点进行一个编号重排
        randomPermutation();
        if (method.equals("LSH")) {
            int k_bins = signatureLength;
            if(n % signatureLength != 0) { k_bins = k_bins + 1; }
            int[] rot_direction = new int[k_bins];
            Random random = new Random();
            for (int i = 0; i < k_bins; i++) {
                if(random.nextBoolean()) { rot_direction[i] = 1; }
                else { rot_direction[i] = -1; }
            }
            // 创建F_OPH数组, 用于存储每个顶点的LSH签名
            OnePermHashSig[] F_OPH_temp = new OnePermHashSig[group_size];
            for (int i = 0; i < group_size; i++) {
                F_OPH_temp[i] = generateSignature(signatureLength, group.get(i), rot_direction);
            }
            // 对分组进行排序
            Integer[] G_temp = new Integer[group_size];
            for (int i = 0; i < group_size; i++) G_temp[i] = i;
            Arrays.sort(G_temp, (o1, o2) -> OnePermHashSig.compare(F_OPH_temp[o1], F_OPH_temp[o2]));
            // 找到第一个组的开始index
            int g_start_temp = 0;
            while (F_OPH_temp[G_temp[g_start_temp]].unassigned())
                g_start_temp++;
            // 计算每个组包含哪些超点
            OnePermHashSig g = new OnePermHashSig(k_bins);
            ArrayList<Integer> Q = new ArrayList<>();
            for (int i = g_start_temp; i < group_size; i++) {
                // 如果遍历到当前顶点的LSH签名与前面的签名不一样，说明前面组的所有顶点都被识别出来，需要增加新的一组了
                if (!F_OPH_temp[G_temp[i]].equals(g)) {
                    // 只返回顶点数量大于1的组
                    if(Q.size() > 1) groups.add(Q);
                    g = F_OPH_temp[G_temp[i]];
                    Q = new ArrayList<>();
                }
                Q.add(group.get(G_temp[i]));
            }
            if(Q.size() > 1) groups.add(Q);
        }
        else if (method.equals("Shingle")) {
            // 初始化F数组, 用于存储每个顶点的shingle值
            int[] F_temp = new int[group_size];
            for (int A = 0; A < group_size; A++)
                F_temp[A] = -1;
            for (int A = 0; A < group_size; A++) {
                // A不是一个超点
                if (I[A] == -1)
                    continue;
                // 将超点A的shingle值先初始化成最大值，然后再逐渐通过超点包含的所有顶点的shingle值逐渐下降
                F_temp[A] = n;
                for (int v = I[A]; ; v = J[v]) {
                    int fv = shingleValue(v);
                    if (F_temp[A] > fv)
                        F_temp[A] = fv;
                    if (J[v] == -1)
                        break;
                }
            }
            // 对分组进行排序
            Integer[] G_temp = new Integer[group_size];
            for (int i = 0; i < group_size; i++) G_temp[i] = i;
            Arrays.sort(G_temp, (o1, o2) -> Integer.compare(F_temp[o1], F_temp[o2]));
            // 找到第一个组的开始index
            int g_start_temp = 0;
            while (F_temp[G_temp[g_start_temp]] == -1)
                g_start_temp++;
            int g = -1;
            ArrayList<Integer> Q = new ArrayList<>();
            for (int i = g_start_temp; i < group_size; i++) {
                if (F_temp[G_temp[i]] != g) {
                    if(Q.size() > 1){
                        groups.add(Q);
                    }
                    g = F_temp[G_temp[i]];
                    Q = new ArrayList<>();
                }
                Q.add(group.get(G_temp[i]));
            }
            if(Q.size() > 1) groups.add(Q);
        }
        return groups;
    }

    /**
     * 顺序切割的分组方案
     * @param iteration 迭代次数
     * @param print_iteration_offset 每隔多少次迭代就进行Encode一次查看压缩率
     */
    public void hierarchicalTest(int iteration, int print_iteration_offset, int maxGroupSize){
        long startTime = System.currentTimeMillis();
        Map<Integer, Map<Integer, Record>> it2Records = new HashMap<>();
        for (int it = 1; it <= iteration; it++) {
            List<List<Integer>> groups = new ArrayList<>();
            logger_.info(String.format("迭代轮数: %d", it));
            // DividePhase =============================================================================================
            {
                logger_.info("开始分组, 采用的是Shingle方法");
                long divideStartTime = System.currentTimeMillis();
                // 首先对顶点进行一个编号重排
                randomPermutation();
                // 初始化F数组, 用于存储每个顶点的shingle值
                int[] F_temp = new int[n];
                for (int A = 0; A < n; A++)
                    F_temp[A] = -1;
                for (int A = 0; A < n; A++) {
                    // A不是一个超点
                    if (I[A] == -1)
                        continue;
                    // 将超点A的shingle值先初始化成最大值，然后再逐渐通过超点包含的所有顶点的shingle值逐渐下降
                    F_temp[A] = n;
                    for (int v = I[A]; ; v = J[v]) {
                        int fv = shingleValue(v);
                        if (F_temp[A] > fv)
                            F_temp[A] = fv;
                        if (J[v] == -1)
                            break;
                    }
                }
                // 对分组进行排序
                Integer[] G_temp = new Integer[n];
                for (int i = 0; i < n; i++) G_temp[i] = i;
                Arrays.sort(G_temp, (o1, o2) -> Integer.compare(F_temp[o1], F_temp[o2]));
                // 找到第一个组的开始index
                int g_start_temp = 0;
                while (F_temp[G_temp[g_start_temp]] == -1)
                    g_start_temp++;

                int g = -1;
                ArrayList<Integer> Q = new ArrayList<>();
                for (int i = g_start_temp; i < n; i++) {
                    if (F_temp[G_temp[i]] != g) {
                        if(Q.size() > 1){
                            if(Q.size() <= maxGroupSize){
                                groups.add(Q);
                            } else {
                                logger_.info(String.format("当前小组顶点数量为%d, 大于阈值%d, 进行层次切割", Q.size(), maxGroupSize));
                                List<List<Integer>> subGroups = hierarchicalSplitGroup(Q, "LSH");
                                for (List<Integer> subGroup : subGroups) {
                                    if(subGroup.size() > 1) groups.add(subGroup);
                                }
                            }
                        }
                        g = F_temp[G_temp[i]];
                        Q = new ArrayList<>();
                    }
                    Q.add(G_temp[i]);
                }
                if(Q.size() > 1) groups.add(Q);
                logger_.info(String.format("分组结束, 花费时间 %5f seconds", (System.currentTimeMillis() - divideStartTime) / 1000.0));
            }
            // End =====================================================================================================

            // MergePhase ==============================================================================================
            {
                double threshold = 1 / ((it + 1) * 1.0);
                logger_.info(String.format("开始合并, 当前顶点的合并阈值:%5f", threshold));
                // 记录mergePhase的运行时间
                long mergeStartTime = System.currentTimeMillis();
                // 记录迭代次数i的所有size>2的组数量
                if(!total_groups_size.containsKey(it))
                    total_groups_size.put(it, new ArrayList<>());
                Map<Integer, Record> recordMap = new HashMap<>();
                // 开始遍历每个组
                for (List<Integer> group : groups) {
                    long start = System.currentTimeMillis();
                    int success = 0;
                    // group_size >= 2时，记录下group_size
                    total_groups_size.get(it).add(group.size());
                    HashMap<Integer, HashMap<Integer, Integer>> hm = createW(group);
                    int initial_size = hm.size();
                    while (hm.size() > 1) {
                        Random rand = new Random();
                        // 从组内随机找到一个超点A
                        int A = rand.nextInt(initial_size);
                        if (hm.get(A) == null) continue;
                        // 变量max记录最大的Jaccard_similarity
                        double max = 0;
                        // 变量idx记录最大Jaccard_similarity的那个顶点
                        int idx = -1;
                        // 遍历组内其他顶点，找到与A的Jaccard Similarity最大的那个顶点
                        for (int j = 0; j < initial_size; j++) {
                            if (hm.get(j) == null)
                                continue;
                            if (j == A) continue;
                            double jaccard_similarity = computeJacSim(hm.get(A), hm.get(j));
                            if (jaccard_similarity > max) {
                                max = jaccard_similarity;
                                idx = j;
                            }
                        }
                        if (idx == -1) {
                            hm.remove(A);
                            continue;
                        }
                        // 这里做了一个交换，目的是把编号较大的顶点合并到编号较小的顶点里面
                        if (group.get(A) > group.get(idx)) {
                            int t = A;
                            A = idx;
                            idx = t;
                        }
                        // 计算两个顶点之间的合并收益
                        double savings = computeSaving(hm.get(A), hm.get(idx), group.get(A), group.get(idx));
                        if (savings >= threshold) {
                            HashMap<Integer, Integer> w_update = updateW(hm.get(A), hm.get(idx));
                            hm.replace(A, w_update);
                            hm.remove(idx);
                            updateSuperNode(group.get(A), group.get(idx));
                            success += 1;
                        } else {
                            hm.remove(A);
                        }
                    }
                    Double time = (System.currentTimeMillis() - start) / 1000.0;
                    if (!recordMap.containsKey(group.size())) {
                        Record r = new Record(group.size(), 1, success, time);
                        recordMap.put(group.size(), r);
                    } else {
                        Record r = recordMap.get(group.size());
                        r.add(1, success, time);
                    }

                }
                it2Records.put(it, recordMap);
                logger_.info(String.format("合并结束, 花费时间 %5f seconds", (System.currentTimeMillis() - mergeStartTime) / 1000.0));
            }
            // End =====================================================================================================
            if (it % print_iteration_offset == 0) {
                logger_.info(String.format("编码结束, 花费时间 %5f seconds", encodePhase_test()));
                evaluatePhase();
                logger_.info(String.format("到这里, 花费时间 %5f seconds", (System.currentTimeMillis() - startTime) / 1000.0));
            }
        }
        for (Integer it : total_groups_size.keySet()) {
            int max = 0;
            for (Integer length : total_groups_size.get(it)) {
                if (max <= length) {
                    max = length;
                }
            }
            logger_.info(String.format("迭代次数 %d: 小组最大顶点数量=%d", it, max));
        }
        logger_.info(String.format("程序运行结束, 总花费时间 %5f seconds", (System.currentTimeMillis() - startTime) / 1000.0));
    }

    /**
     * @param iteration              迭代次数
     * @param print_iteration_offset 每执行多少次迭代就进行一次 encode 和 evaluate 进行结果输出
     */
    @Override
    public void run(int iteration, int print_iteration_offset) {
//        originTest(iteration, print_iteration_offset);
//        sequentialTest(iteration, print_iteration_offset, 2000);
        hierarchicalTest(iteration, print_iteration_offset, 500);
//        long starTime = System.currentTimeMillis();
//        for (int it = 1; it <= iteration; it++) {
//            logger_.info(String.format("迭代轮数: %d", it));
//            double threshold = 1 / ((it + 1) * 1.0);
//            logger_.info(String.format("分组结束, 花费时间 %5f seconds", dividePhase()));
//            System.out.println(String.format("@Time: %5f seconds", mergePhase(threshold)));
//            logger_.info(String.format("合并结束, 花费时间 %5f seconds", mergePhase_test(it, threshold)));
//            if (it % print_iteration_offset == 0) {
//                logger_.info(String.format("编码结束, 花费时间 %5f seconds", encodePhase_test()));
//                evaluatePhase();
//                logger_.info(String.format("到这里, 花费时间 %5f seconds", (System.currentTimeMillis() - starTime) / 1000.0));
//            }
//        }
//
//        for (Integer it : total_groups_size.keySet()) {
//            int max = 0;
//            for (Integer length : total_groups_size.get(it)) {
//                if (max <= length) {
//                    max = length;
//                }
//            }
//            logger_.info(String.format("迭代次数 %d: 小组最大顶点数量=%d", it, max));
//        }
//        logger_.info(String.format("程序运行结束, 总花费时间 %5f seconds", (System.currentTimeMillis() - starTime) / 1000.0));
    }
}
