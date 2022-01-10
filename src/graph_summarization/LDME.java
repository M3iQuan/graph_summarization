package graph_summarization;

import org.javatuples.Pair;

import java.io.*;
import java.util.*;
import java.util.concurrent.ThreadLocalRandom;
import java.util.zip.Adler32;

public class LDME extends Summary{

    // 用于对顶点的重新编号
    int[] h;
    // 哈希签名的长度
    int signatureLength;
    // 哈希签名数组，在顶点分组时使用
    OnePermHashSig[] F_OPH;
    Integer[] G;
    // 找到G数组里的第一个组的第一个顶点的index
    int g_start;
    // 分组后的组别数量
    int num_groups;
    // 用于分组的两个数组，分别记录每个组的hash值和第几组
    OnePermHashSig[] group_prop_0;
    int[] group_prop_1;

    // 用于计算合并的成功率
    ArrayList<Double> merged_success;
    // 用于统计每次迭代时分组的顶点数量
    HashMap<Integer, ArrayList<Integer>> total_groups_size;


    /**
     * 构造函数，用于初始化一些共同的结构
     *
     * @param basename 数据集的基本名字
     * @throws Exception 抛出异常
     */
    public LDME(String basename, int signatureLength) throws Exception {
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
     * 计算一个组包含多少的顶点，即组的大小
     *
     * @param arr    排序后各组的哈希值，按照哈希值升序排序
     * @param key    要计算大小的那组的哈希值
     * @param st_pos 要计算大小的那组的开始index
     * @return
     */
    private int groupLength(OnePermHashSig[] arr, OnePermHashSig key, int st_pos) {
        int counter = 1;
        while (arr[st_pos++].equals(key) && st_pos < arr.length - 1) counter++;
        return counter;
    }

    /**
     * 分组阶段，LDME算法通过哈希值对顶点进行划分
     * 分组完成后，可以通过遍历 G 得到每个组，相同组的G[i]值相等
     */
    @Override
    public double dividePhase(){
        long startTime = System.currentTimeMillis();
//        System.out.println("# Divide Phase");
        logger_.info("开始分组, 采用的是Local Sensitive Hash方法");
        int k_bins = signatureLength;
        int bin_size = n / k_bins;
        if (n % k_bins != 0) { k_bins = k_bins + 1; }

        // 首先生成长度为k_bins的一个数组用于辅助计算hash签名值
        int[] rot_direction = new int[k_bins];
        Random random = new Random();
        for (int i = 0; i < k_bins; i++) {
            if (random.nextBoolean()) { rot_direction[i] = 1; }
            else { rot_direction[i] = -1; }
        }

        // 接着对顶点进行一个编号重排
        randomPermutation();

        // 初始化F_OPH数组, 用于存储每个顶点的哈希值
        F_OPH = new OnePermHashSig[n];
        for(int A=0; A<n ; A++)
            F_OPH[A] = new OnePermHashSig(k_bins);

        for (int A = 0; A < n; A++) {
            // A不是一个超点
            if (I[A] == -1) continue;
            // 遍历超点的每个顶点v
            for (int v = I[A]; ; v=J[v]) {
                int[] neighbours = neighbors_[v];
//                int[] neighbours = Gr.successorArray(v);
                // 遍历顶点v的每个邻居j
                for (int j = 0; j < neighbours.length; j++) {
                    // 得到邻居j重排后的id
                    int permuted_h = h[neighbours[j]];
                    // 以及在第几个bin
                    int permuted_bin = permuted_h / bin_size;
                    // 如果b_i==-1(即还没设置)或者比b_i还小的索引 (permuted_h%bin_size) 存在顶点，则重新设置b_i
                    if (F_OPH[A].sig[permuted_bin] == -1 || permuted_h % bin_size < F_OPH[A].sig[permuted_bin]) {
                        F_OPH[A].sig[permuted_bin] = permuted_h % bin_size;
                    }
                }
                // 遍历完超点A的所有顶点就直接退出
                if(J[v]==-1)
                    break;
            }

            // rotation
            for (int A_bin = 0; A_bin < k_bins; A_bin++) {
                int direction = rot_direction[A_bin];
                // 如果b_i还没设置, 则利用rot_direction来确定往哪边采样
                if (F_OPH[A].sig[A_bin] == -1) {
                    int i = (A_bin + direction) % k_bins;
                    if (i < 0) { i += k_bins; }
                    int counter = 0;
                    while (F_OPH[A].sig[i] == -1 && counter < k_bins) {
                        i = (i + direction) % k_bins;
                        if (i < 0) { i += k_bins; }
                        counter++;
                    }
                    F_OPH[A].sig[A_bin] = F_OPH[A].sig[i];
                }
            }

        }

        // 对分组进行排序
        G = new Integer[n];
        for (int i = 0; i < n; i++) G[i] = i;
        Arrays.sort(G, (o1, o2) -> OnePermHashSig.compare(F_OPH[o1], F_OPH[o2]));

        // 找到第一个组的开始index
        g_start = 0;
        while (F_OPH[G[g_start]].unassigned())
            g_start++;

        // 计算组的数量
        num_groups = 0;
        OnePermHashSig g = new OnePermHashSig(k_bins);
        for (int i = g_start; i < n; i++) {
            if (!F_OPH[G[i]].equals(g)) {
                num_groups++;
                g = F_OPH[G[i]];
            }
        }

        // 当F[G[i]]>F[G[i-1]] 形成新的组
        group_prop_0 = new OnePermHashSig[num_groups];
        group_prop_1 = new int[num_groups];
        g = new OnePermHashSig(k_bins + 1);
        int counter = 0;
        for (int i = g_start; i < n; i++) {
            if (!F_OPH[G[i]].equals(g)) {
                group_prop_0[counter] = F_OPH[G[i]];
                group_prop_1[counter] = i;
                counter++;
                g = F_OPH[G[i]];
            }
        }

//        System.out.println("g_start:" + g_start + ", num_groups:" + num_groups);
//        logger_.info(String.format("分组结束, num_groups:%d", num_groups));
        return (System.currentTimeMillis() - startTime) / 1000.0;
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
     * 传入参数k，对所有顶点进行分组
     * @param k LSH签名的长度
     * @return 返回生成的组，每一个ArrayList都是一组，组内元素是超点的编号
     */
    public List<List<Integer>> dividePhase_new(int k){
        System.out.println("# Divide Phase");

        int k_bins = k;
        if (n % k_bins != 0) { k_bins = k_bins + 1; }

        // 首先生成长度为k_bins的一个数组用于辅助计算hash签名值
        int[] rot_direction = new int[k_bins];
        Random random = new Random();
        for (int i = 0; i < k_bins; i++) {
            if (random.nextBoolean()) { rot_direction[i] = 1; }
            else { rot_direction[i] = -1; }
        }

        // 接着对顶点进行一个编号重排
        randomPermutation();

        // 创建F_OPH数组, 用于存储每个顶点的LSH签名
        OnePermHashSig[] F_OPH_temp = new OnePermHashSig[n];
        for (int A = 0; A < n; A++) {
            F_OPH_temp[A] = generateSignature(k, A, rot_direction);
        }

        // 对分组进行排序
        Integer[] G_temp = new Integer[n];
        for (int i = 0; i < n; i++) G_temp[i] = i;
        Arrays.sort(G_temp, (o1, o2) -> OnePermHashSig.compare(F_OPH_temp[o1], F_OPH_temp[o2]));

        // 找到第一个组的开始index
        int g_start_temp = 0;
        while (F_OPH_temp[G_temp[g_start_temp]].unassigned())
            g_start_temp++;

        // 计算每个组包含哪些超点
        OnePermHashSig g = new OnePermHashSig(k_bins);
        List<List<Integer>> groups = new ArrayList<>();
        ArrayList<Integer> Q = new ArrayList<>();
        for (int i = g_start_temp; i < n; i++) {
            // 如果遍历到当前顶点的LSH签名与前面的签名不一样，说明前面组的所有顶点都被识别出来，需要增加新的一组了
            if (!F_OPH_temp[G_temp[i]].equals(g)) {
                // 只返回顶点数量大于1的组
                if(Q.size() > 1) groups.add(Q);
                g = F_OPH_temp[G_temp[i]];
                Q = new ArrayList<>();
            }
            Q.add(G_temp[i]);
        }

        return groups;
    }

    /**
     * 对于某个小组内的顶点进行分组
     * @param k LSH签名的长度
     * @param group 需要再进行划分的小组
     * @return 返回生成的组，每一个ArrayList都是一组，组内元素是超点的编号
     */
    public List<List<Integer>> hierarchicalDivide(int k, List<Integer> group) {
//        System.out.println("# Hierarchical Divide Phase");

        int group_size = group.size();

        int k_bins = k;
        if (n % k_bins != 0) { k_bins = k_bins + 1; }

        // 首先生成长度为k_bins的一个数组用于辅助计算hash签名值
        int[] rot_direction = new int[k_bins];
        Random random = new Random();
        for (int i = 0; i < k_bins; i++) {
            if (random.nextBoolean()) { rot_direction[i] = 1; }
            else { rot_direction[i] = -1; }
        }

        // 接着对顶点进行一个编号重排
        randomPermutation();

        // 创建F_OPH数组, 用于存储每个顶点的LSH签名
        OnePermHashSig[] F_OPH_temp = new OnePermHashSig[group_size];
        for (int i = 0; i < group_size; i++) {
            F_OPH_temp[i] = generateSignature(k, group.get(i), rot_direction);
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
        List<List<Integer>> result = new ArrayList<>();
        ArrayList<Integer> Q = new ArrayList<>();
        for (int i = g_start_temp; i < group_size; i++) {
            // 如果遍历到当前顶点的LSH签名与前面的签名不一样，说明前面组的所有顶点都被识别出来，需要增加新的一组了
            if (!F_OPH_temp[G_temp[i]].equals(g)) {
                // 只返回顶点数量大于1的组
                if(Q.size() > 1) result.add(Q);
                g = F_OPH_temp[G_temp[i]];
                Q = new ArrayList<>();
            }
            Q.add(group.get(G_temp[i]));
        }
        if(Q.size() > 1) result.add(Q);

        return result;
    }


    public List<List<Integer>> hierarchicalDivide(List<Integer> group){
        HashCode[] hashCode = new HashCode[group.size()];
        Integer[] G_temp = new Integer[group.size()];
        for (int i = 0; i < group.size(); i++) {
            G_temp[i] = i;
            hashCode[i] = generateHashCode_2(group.get(i));
        }
        Arrays.sort(G_temp, (o1, o2) -> HashCode.compare(hashCode[o1], hashCode[o2]));

        int g_start_temp = 0;
        List<List<Integer>> result = new ArrayList<>();
        int k_bins = 1;
        HashCode g_code = new HashCode(k_bins);
        ArrayList<Integer> Q = new ArrayList<>();
        for (int i = g_start_temp; i < group.size(); i++) {
            // 如果遍历到当前顶点的LSH签名与前面的签名不一样，说明前面组的所有顶点都被识别出来，需要增加新的一组了
            if (!hashCode[G_temp[i]].equals(g_code)) {
                // 只返回顶点数量大于1的组
                if (Q.size() > 1) result.add(Q);
                g_code = hashCode[G_temp[i]];
                Q = new ArrayList<>();
            }
            Q.add(group.get(G_temp[i]));
        }
        if(Q.size() > 1) result.add(Q);

        return result;
    }

    public double divideInGroups(int[] Q, int group_size){
        long startTime = System.currentTimeMillis();
        System.out.println("# Divide Phase");

        int k_bins = signatureLength;
        int bin_size = n / k_bins;
        if (n % k_bins != 0) { k_bins = k_bins + 1; }

        // 首先生成长度为k_bins的一个数组用于辅助计算hash签名值
        int[] rot_direction = new int[k_bins];
        Random random = new Random();
        for (int i = 0; i < k_bins; i++) {
            if (random.nextBoolean()) { rot_direction[i] = 1; }
            else { rot_direction[i] = -1; }
        }

        // 接着对顶点进行一个编号重排
        randomPermutation();

        // 初始化F_OPH_temp数组, 用于存储每个顶点的哈希值
        OnePermHashSig[] F_OPH_temp = new OnePermHashSig[n];
        for(int A=0; A<n ; A++)
            F_OPH_temp[A] = new OnePermHashSig(k_bins);

        TreeSet<Integer> treeSet = new TreeSet<>();
        for (int i = 0; i < Q.length; i++) {
            treeSet.add(Q[i]);
        }

        for (int A = 0; A < n; A++) {
            // A不是一个超点
            if (I[A] == -1) continue;
            if(!treeSet.contains(I[A])) continue;
            // 遍历超点的每个顶点v
            for (int v = I[A]; ; v=J[v]) {
                int[] neighbours = neighbors_[v];
//                int[] neighbours = Gr.successorArray(v);
                // 遍历顶点v的每个邻居j
                for (int j = 0; j < neighbours.length; j++) {
                    // 得到邻居j重排后的id
                    int permuted_h = h[neighbours[j]];
                    // 以及在第几个bin
                    int permuted_bin = permuted_h / bin_size;
                    // 如果b_i==-1(即还没设置)或者比b_i还小的索引 (permuted_h%bin_size) 存在顶点，则重新设置b_i
                    if (F_OPH[A].sig[permuted_bin] == -1 || permuted_h % bin_size < F_OPH[A].sig[permuted_bin]) {
                        F_OPH[A].sig[permuted_bin] = permuted_h % bin_size;
                    }
                }
                // 遍历完超点A的所有顶点就直接退出
                if(J[v]==-1)
                    break;
            }

            // rotation
            for (int A_bin = 0; A_bin < k_bins; A_bin++) {
                int direction = rot_direction[A_bin];
                // 如果b_i还没设置, 则利用rot_direction来确定往哪边采样
                if (F_OPH[A].sig[A_bin] == -1) {
                    int i = (A_bin + direction) % k_bins;
                    if (i < 0) { i += k_bins; }
                    int counter = 0;
                    while (F_OPH[A].sig[i] == -1 && counter < k_bins) {
                        i = (i + direction) % k_bins;
                        if (i < 0) { i += k_bins; }
                        counter++;
                    }
                    F_OPH[A].sig[A_bin] = F_OPH[A].sig[i];
                }
            }

        }

        // 对分组进行排序
        Integer[] G_temp = new Integer[n];
        for (int i = 0; i < n; i++) G[i] = i;
        Arrays.sort(G, (o1, o2) -> OnePermHashSig.compare(F_OPH[o1], F_OPH[o2]));

        // 找到第一个组的开始index
        g_start = 0;
        while (F_OPH[G[g_start]].unassigned())
            g_start++;

        // 计算组的数量
        num_groups = 0;
        OnePermHashSig g = new OnePermHashSig(k_bins);
        for (int i = g_start; i < n; i++) {
            if (!F_OPH[G[i]].equals(g)) {
                num_groups++;
                g = F_OPH[G[i]];
            }
        }

        // 当F[G[i]]>F[G[i-1]] 形成新的组
        group_prop_0 = new OnePermHashSig[num_groups];
        group_prop_1 = new int[num_groups];
        g = new OnePermHashSig(k_bins + 1);
        int counter = 0;
        for (int i = g_start; i < n; i++) {
            if (!F_OPH[G[i]].equals(g)) {
                group_prop_0[counter] = F_OPH[G[i]];
                group_prop_1[counter] = i;
                counter++;
                g = F_OPH[G[i]];
            }
        }

        return (System.currentTimeMillis() - startTime) / 1000.0;
    }

    /**
     * 合并阶段，LDME算法在每个小组内的顶点采用Random方式进行合并
     */
    @Override
    public double mergePhase(double threshold){
        System.out.println("# Merge Phase");
        System.out.printf("Threshold=%5f%n", threshold);
        long startTime = System.currentTimeMillis();
        int idx = 0;
        OnePermHashSig[] temp = new OnePermHashSig[n];
        for (int i = 0; i < n; i++) temp[i] = F_OPH[G[i]];

        // 开始遍历每个组
        for (int i = 0; i < num_groups; i++) {
            int st_position = group_prop_1[i];
            int group_size = groupLength(temp, group_prop_0[i], st_position) - 1;
            // 如果一个组只有一个顶点则直接跳过该组
            if (group_size < 2) continue;

            // Q是当前要进行合并的的组
            int[] Q = new int[group_size];
            int counter = 0;
            // 找到在Q组里面的所有超点
            for (int j= st_position; j < (st_position + group_size); j++){
                Q[counter++] = G[j];
            }

            HashMap<Integer, HashMap<Integer,Integer>> hm = createW(Q, group_size);
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
    public double mergePhase_test(double threshold){
        System.out.println("# Merge Phase");
        System.out.printf("Threshold=%5f%n", threshold);
        long startTime = System.currentTimeMillis();
        int idx = 0;
        OnePermHashSig[] temp = new OnePermHashSig[n];
        for (int i = 0; i < n; i++) temp[i] = F_OPH[G[i]];

        HashMap<Integer, ArrayList<Integer>> merged_ = new HashMap<>();

        // 开始遍历每个组
        for (int i = 0; i < num_groups; i++) {
            int st_position = group_prop_1[i];
            int group_size = groupLength(temp, group_prop_0[i], st_position) - 1;
            // 如果一个组只有一个顶点则直接跳过该组
            if (group_size < 2) continue;

            // 能成功合并的对数
            int merged_pairs = 0;
            if (!merged_.containsKey(group_size)) {
                merged_.put(group_size, new ArrayList<>());
            }

            // Q是当前要进行合并的的组
            int[] Q = new int[group_size];
            int counter = 0;
            // 找到在Q组里面的所有超点
            for (int j= st_position; j < (st_position + group_size); j++){
                Q[counter++] = G[j];
            }

            HashMap<Integer, HashMap<Integer,Integer>> hm = createW(Q, group_size);
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
     * 测试函数，用于验证一下想法：
     *     1. 每次迭代时小组顶点数量
     *
     */
    public double mergePhase_test2(int iteration, double threshold){
        logger_.info(String.format("开始合并, 当前顶点的合并阈值:%5f", threshold));
//        System.out.println("# Merge Phase");
//        System.out.printf("Threshold=%5f%n", threshold);
        // 记录mergePhase的运行时间
        long startTime = System.currentTimeMillis();
        int idx = 0;
        OnePermHashSig[] temp = new OnePermHashSig[n];
        for (int i = 0; i < n; i++) temp[i] = F_OPH[G[i]];

        // 记录下 <group_size, merged_pairs>， 用来计算每个组内的成功合并率 merged_pairs/group_size
        HashMap<Integer, ArrayList<Integer>> merged_ = new HashMap<>();

        // 记录迭代次数i的所有size>2的组数量
        if(!total_groups_size.containsKey(iteration))
            total_groups_size.put(iteration, new ArrayList<>());

        // 开始遍历每个组
        for (int i = 0; i < num_groups; i++) {
            int st_position = group_prop_1[i];
            int group_size = groupLength(temp, group_prop_0[i], st_position) - 1;

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
            for (int j= st_position; j < (st_position + group_size); j++){
                Q[counter++] = G[j];
            }

            HashMap<Integer, HashMap<Integer,Integer>> hm = createW(Q, group_size);
            int initial_size = hm.size();
            while (hm.size() > 1) {
                Random rand = new Random();
                // 从组内随机找到一个超点A
                int A = rand.nextInt(initial_size);
                if (hm.get(A) == null)
                    continue;

                // 变量max记录最大的jaccard_similarity
                double max = 0;
                // 变量idx记录最大jaccard_similarity的那个顶点
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
     * 通过dividePhase_new()函数进行分组后，产生的groups传入这个函数进行按小组合并
     * @param threshold 合并的阈值
     * @param groups 所有的小组
     * @return 合并所需要的所有时间
     */
    public double mergePhase_new(double threshold, List<List<Integer>> groups){
        System.out.println("# Merge Phase");
        System.out.printf("Threshold=%5f%n", threshold);
        // 记录mergePhase的运行时间
        long startTime = System.currentTimeMillis();

        // 开始遍历每个组
        for (List<Integer> Q : groups) {
            mergeInGroup(Q, threshold);
        }

        return (System.currentTimeMillis() - startTime) / 1000.0;
    }


    /**
     * 对组内数量大于阈值max_group_size的小组进行顺序切分后再按小组合并
     * @param threshold 合并的阈值
     * @param groups 所有的小组
     * @return 合并所需要的时间
     */
    public double mergePhaseSequentialSplit(double threshold, List<List<Integer>> groups) {
        System.out.println("# Merge Phase");
        System.out.printf("Threshold=%5f%n", threshold);
        int max_group_size = 1000;
        // 记录mergePhase的运行时间
        long startTime = System.currentTimeMillis();

        // 开始遍历每个组
        for (List<Integer> Q : groups) {
            int group_size = Q.size();
            // 如果组内顶点数量比较少
            if (group_size <= max_group_size) {
                mergeInGroup(Q, threshold);
                continue;
            }
            List<List<Integer>> subGroups = sequentialSplitGroup(Q, max_group_size);
            for (List<Integer> subGroup : subGroups) {
                mergeInGroup(subGroup, threshold);
            }
        }
        return (System.currentTimeMillis() - startTime) / 1000.0;
    }


    public double mergePhaseHierarchicalDivide(double threshold, List<List<Integer>> groups){
        System.out.println("# Merge Phase");
        System.out.printf("Threshold=%5f%n", threshold);
        int max_group_size = 1000;
        // 记录mergePhase的运行时间
        long startTime = System.currentTimeMillis();

        // 开始遍历每个组
        for (List<Integer> Q : groups) {
            int group_size = Q.size();
            // 如果组内顶点数量比较少
            if (group_size <= max_group_size) {
                mergeInGroup(Q, threshold);
                continue;
            }
            List<List<Integer>> subGroups = hierarchicalDivide(2*signatureLength, Q);
            for (List<Integer> subGroup : subGroups) {
                mergeInGroup(subGroup, threshold);
            }
        }
        return (System.currentTimeMillis() - startTime) / 1000.0;
    }


    /**
     * 直接对组Q内的超点按Random算法进行合并
     * @param Q 小组Q通过int[]输入
     * @param group_size 小组Q的顶点数量
     * @param threshold 合并阈值
     */
    public void mergeInGroup(int[] Q, int group_size, double threshold) {
        int idx = 0;
        HashMap<Integer, HashMap<Integer,Integer>> hm = createW(Q, group_size);
        int initial_size = hm.size();
        while (hm.size() > 1) {
            Random rand = new Random();
            // 从组内随机找到一个超点A
            int A = rand.nextInt(initial_size);
            if (hm.get(A) == null)
                continue;

            // 变量max记录最大的jaccard_similarity
            double max = 0;
            // 变量idx记录最大jaccard_similarity的那个顶点
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

    /**
     * 直接对组Q内的超点按Random算法进行合并
     * @param Q 小组Q通过List<>输入
     * @param threshold 合并阈值
     */
    public void mergeInGroup(List<Integer> Q, double threshold){
        HashMap<Integer, HashMap<Integer,Integer>> hm = createW(Q);
        int initial_size = hm.size();
        while (hm.size() > 1) {
            Random rand = new Random();
            // 从组内随机找到一个超点A
            int A = rand.nextInt(initial_size);
            if (hm.get(A) == null)
                continue;

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
            if (Q.get(A) > Q.get(idx)) {
                int t = A;
                A = idx;
                idx = t;
            }

            // 计算两个顶点之间的合并收益
            double savings = computeSaving(hm.get(A), hm.get(idx), Q.get(A), Q.get(idx));
            if (savings >= threshold) {
                HashMap<Integer, Integer> w_update = updateW(hm.get(A), hm.get(idx));
                hm.replace(A, w_update);
                hm.remove(idx);
                updateSuperNode(Q.get(A), Q.get(idx));
            } else {
                hm.remove(A);
            }
        }

    }

    /**
     * 对组内顶点数量大于maxGroupSize的小组进行顺序划分
     * @param group 组内数量大于maxGroupSize的小组
     * @param maxGroupSize 最大组内数量阈值
     * @return 按阈值切分后的小组
     */
    public List<List<Integer>> sequentialSplitGroup(List<Integer> group, int maxGroupSize) {
        List<List<Integer>> groups = new ArrayList<>();
        int group_size = group.size();
        int num = group_size / maxGroupSize;
        if (group_size % maxGroupSize != 0) num += 1;
        for (int j = 0; j < num; j++) {
            int start = j * maxGroupSize;
            int end = Math.min((j + 1) * maxGroupSize, group_size);
            // Q是当前要进行合并的的组
            List<Integer> subGroup = new ArrayList<>();
            for (int k = start; k < end; k++) {
                subGroup.add(group.get(k));
            }
            groups.add(subGroup);
        }
        return groups;
    }

    /**
     *
     * @param iteration
     * @param threshold
     * @return
     */
    public double multiDivide(int iteration, double threshold){
        System.out.println("# Multi Divide ");
        System.out.printf("Threshold=%5f%n", threshold);
        // 记录mergePhase的运行时间
        long startTime = System.currentTimeMillis();
        int idx = 0;

        // 第一次分组，记住所有顶点的分组id
        System.out.printf("First Divide: %5f seconds%n", dividePhase());
        OnePermHashSig[] temp1 = new OnePermHashSig[n];
        for (int i = 0; i < n; i++) temp1[i] = F_OPH[G[i]];
        Map<Integer, Integer> node_to_group = new HashMap<>();
        for (int i = 0; i < num_groups; i++) {
            int st_position = group_prop_1[i];
            int group_size = groupLength(temp1, group_prop_0[i], st_position) - 1;
            TreeSet<Integer> Q = new TreeSet<>();
            // 找到在Q组里面的所有超点
            for (int j= st_position; j < (st_position + group_size); j++){
                node_to_group.put(G[j], i);
            }
        }

        // 第二次分组，加多一个限制，只有组内超点的第一次分组在同一个组才能真正在同一个组
        System.out.printf("Second Divide: %5f seconds%n", dividePhase());
        OnePermHashSig[] temp2 = new OnePermHashSig[n];
        int origin_num = 0;
        int update_num = 0;
        int update_num_minus_1 = 0;
        int origin_max = 0;
        int update_max = 0;

        for (int i = 0; i < n; i++) temp2[i] = F_OPH[G[i]];
        for (int i = 0; i < num_groups; i++) {
            int st_position = group_prop_1[i];
            int group_size = groupLength(temp2, group_prop_0[i], st_position) - 1;
            if(group_size < 2) continue;
            TreeSet<Integer> Q = new TreeSet<>();
            // 找到在Q组里面的所有超点
            if(group_size > 500)
                System.out.println();
            if(group_size > 500) System.out.print("Group " + i + " contains " + group_size + " nodes:");
            for (int j= st_position; j < (st_position + group_size); j++){
                Q.add(G[j]);
                if(group_size > 500)
                    System.out.print(G[j] + " ");
            }
            if(group_size > 500)
                System.out.println();
            origin_num += 1;
            if(group_size > origin_max) origin_max = group_size;

            HashMap<Integer, TreeSet<Integer>> group_to_node = new HashMap<>();
            for (Integer A : Q) {
                if (!group_to_node.containsKey(node_to_group.get(A))) {
                    group_to_node.put(node_to_group.get(A), new TreeSet<>());
                }
                group_to_node.get(node_to_group.get(A)).add(A);
            }

            for (Integer j : group_to_node.keySet()) {
                update_num += 1;
                if(group_to_node.get(j).size() > 2) update_num_minus_1 += 1;
                if(group_to_node.get(j).size() > update_max) update_max = group_to_node.get(j).size();
                if(group_size > 500){
                    System.out.print(group_to_node.get(j).size() + "(");
                    for(Integer A : group_to_node.get(j)){
                        System.out.print(A+",");
                    }
                    System.out.print(") ");
                }


            }
            if(group_size > 500)
                System.out.println();
        }
        System.out.println("Origin groups num:" + origin_num + ", maximum_group_size:" + origin_max);
        System.out.println("Update groups num:" + update_num + ", maximum_group_size:" + update_max);
        System.out.println("Update groups bigger than 1 num:" + update_num_minus_1 + ", maximum_group_size:" + update_max);
//        File writeFile = new File("/home/hallo/下载/GraphSummarization/result/20211102/multiDivide_5.txt");
//        try {
//            if(!writeFile.exists()) writeFile.createNewFile();
//            BufferedWriter out = new BufferedWriter(new FileWriter(writeFile));
//            int[] count= new int[n];
//            for(int i=0; i <n; i++) count[i]=0;
//
//            for (Integer A : divide_1.keySet()) {
//                count[A] += 1;
//                out.write(A + " ===============================\n");
//                out.write("1: ");
//                for(Integer B : divide_1.get(A)){
//                    out.write(B + " ");
//                }
//                out.write("\n");
//                if(divide_2.containsKey(A)) {
//                    out.write("2: ");
//                    for (Integer B : divide_2.get(A)) {
//                        out.write(B + " ");
//                    }
//                    out.write("\n");
//                }
//                out.write("=====================================\n");
//                out.flush();
//            }
//
//            for (Integer A : divide_2.keySet()) {
//                if(count[A] != 0) continue;
//                out.write(A + " ===============================\n");
//                if(divide_1.containsKey(A)){
//                    out.write("1: ");
//                    for(Integer B : divide_1.get(A)){
//                        out.write(B + " ");
//                    }
//                    out.write("\n");
//                }
//                out.write("2: ");
//                for(Integer B : divide_2.get(A)){
//                    out.write(B + " ");
//                }
//                out.write("\n");
//                out.write("=====================================\n");
//                out.flush();
//            }
//
//        } catch (Exception e) {
//            System.out.println(e);
//        }

        return (System.currentTimeMillis() - startTime) / 1000.0;
    }

    /**
     * 输出结果
     *
     */
    public void outputResult(Map<Integer, Map<Integer, Record>> it2Records, String method, int counter) {
        String filename = "/home/hallo/下载/GraphSummarization/result/20211126/" + method + "_" + signatureLength + "_" + counter +".txt";
        File writeFile = new File(filename);
        // 遍历每一次迭代
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(writeFile));
            for (Integer it : new TreeSet<>(it2Records.keySet())) {
                Map<Integer, Record> recordMap = it2Records.get(it);
                Set<Integer> order = new TreeSet<>(recordMap.keySet());
                for (Integer group_size : order) {
                    Record r = recordMap.get(group_size);
                    out.write(r.toString() + " ");
                    out.flush();
                }
                out.write("\n");
                out.flush();
            }
            out.close();
        } catch (Exception e) {
            System.out.println(e.getMessage());
        }
    }

    /**
     * 把迭代过程中得到的每个小组合并的成功次数输出出来， 每一行记录的是 group_size1,num1,success1 group_size2,num2,success2 ...
     * @param it2MergedSuccess 记录了每次迭代的小组合并的成功次数
     * @param method 进行实验的方法： origin, sequential 和 hierarchical.txt
     * @param counter 第几次实验
     */
    public void analysisMerged(Map<Integer, Map<Integer, Pair<Integer, Integer>>> it2MergedSuccess, String method, int counter){
        String filename = "/home/hallo/下载/GraphSummarization/result/20211127/" + method + "_mergedSuccess_" + counter +".txt";
        File writeFile = new File(filename);
        // 遍历每一次迭代
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(writeFile));
            for (Integer it : it2MergedSuccess.keySet()) {
                Map<Integer, Pair<Integer, Integer>> mergedSuccess = it2MergedSuccess.get(it);
                Set<Integer> order = new TreeSet<>(mergedSuccess.keySet());
                for (Integer group_size : order) {
                    Integer num = mergedSuccess.get(group_size).getValue0();
                    Integer success = mergedSuccess.get(group_size).getValue1();
                    out.write(group_size + "," + num + "," + success + " ");
                    out.flush();
                }
                out.write("\n");
                out.flush();
            }
            out.close();
        } catch (Exception e) {
            System.out.println(e.getMessage());
        }
    }

    /**
     * 记录origin方法中小组的大小分布和合并成功次数
     * @param iteration 迭代次数
     * @param print_iteration_offset 每隔多少次迭代就进行Encode一次查看压缩率
     */
    public void originTest(int iteration, int print_iteration_offset) {
        System.out.println("----------------------------------- ORIGIN ALGORITHM ---------------------------------------");
        Map<Integer, Map<Integer, Record>> it2Records = new HashMap<>();
        for (int it = 1; it <= iteration; it++) {
            List<List<Integer>> groups = new ArrayList<>();
            System.out.println("\n------------------------- ITERATION " + it);
            // DividePhase =============================================================================================
            {
                System.out.println("# Divide Phase");
                long divideStartTime = System.currentTimeMillis();
                int k_bins = signatureLength;
                if (n % k_bins != 0) {
                    k_bins = k_bins + 1;
                }
                // 首先生成长度为k_bins的一个数组用于辅助计算hash签名值
                int[] rot_direction = new int[k_bins];
                Random random = new Random();
                for (int i = 0; i < k_bins; i++) {
                    if (random.nextBoolean()) {
                        rot_direction[i] = 1;
                    } else {
                        rot_direction[i] = -1;
                    }
                }
                // 接着对顶点进行一个编号重排
                randomPermutation();
                // 创建F_OPH数组, 用于存储每个顶点的LSH签名
                OnePermHashSig[] F_OPH_temp = new OnePermHashSig[n];
                for (int A = 0; A < n; A++) {
                    F_OPH_temp[A] = generateSignature(signatureLength, A, rot_direction);
                }
                // 对分组进行排序
                Integer[] G_temp = new Integer[n];
                for (int i = 0; i < n; i++) G_temp[i] = i;
                Arrays.sort(G_temp, (o1, o2) -> OnePermHashSig.compare(F_OPH_temp[o1], F_OPH_temp[o2]));
                // 找到第一个组的开始index
                int g_start_temp = 0;
                while (F_OPH_temp[G_temp[g_start_temp]].unassigned())
                    g_start_temp++;
                // 计算每个组包含哪些超点
                OnePermHashSig g = new OnePermHashSig(k_bins);
                ArrayList<Integer> Q = new ArrayList<>();
                for (int i = g_start_temp; i < n; i++) {
                    // 如果遍历到当前顶点的LSH签名与前面的签名不一样，说明前面组的所有顶点都被识别出来，需要增加新的一组了
                    if (!F_OPH_temp[G_temp[i]].equals(g)) {
                        // 只返回顶点数量大于1的组
                        if (Q.size() > 1) groups.add(Q);
                        g = F_OPH_temp[G_temp[i]];
                        Q = new ArrayList<>();
                    }
                    Q.add(G_temp[i]);
                }
                if(Q.size() > 1) groups.add(Q);
                System.out.printf("@Time: %5f seconds%n", (System.currentTimeMillis() - divideStartTime) / 1000.0);
            }
            // End =====================================================================================================

            // MergePhase ==============================================================================================
            {
                double threshold = 1 / ((it + 1) * 1.0);
                System.out.println("# Merge Phase");
                System.out.printf("Threshold=%5f%n", threshold);
                // 记录mergePhase的运行时间
                long mergeStartTime = System.currentTimeMillis();
                // 开始遍历每个组
                Map<Integer, Record> recordMap = new HashMap<>();
                for (List<Integer> group : groups) {
                    long start = System.currentTimeMillis();
                    int success = 0;
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
                System.out.printf("@Time: %5f seconds%n", (System.currentTimeMillis() - mergeStartTime) / 1000.0);
            }
            // End =====================================================================================================
            if (it % print_iteration_offset == 0) {
                System.out.printf("@Time: %5f seconds%n", encodePhase_new());
                evaluatePhase();
            }
        }
        outputResult(it2Records, "origin", 1);
    }

    /**
     * 对顶点生成一个标识，<neighbors_num, edges_num> 用来进行细分
     * @param super_node_id
     * @return
     */
    public HashCode generateHashCode_2(int super_node_id) {
        int k_bins = 1;
        HashCode hashCode = new HashCode(k_bins);
        Map<Integer, Integer> singleW = createSingleW(super_node_id);
        for (Integer u : singleW.keySet()) {
            hashCode.code[0] += 1;  // 邻居的数量
//            hashCode.code[1] += singleW.get(u);  // 与邻居的总边数， 除以超点大小就得到了连接权重
        }
        return hashCode;
    }

    public HashCodeAdjust generateHashCode_3(int super_node_id){
        int k_bins = 2;
        int superNodeSize = superNodeLength(super_node_id);
        HashCodeAdjust hashCode = new HashCodeAdjust(k_bins);
        Map<Integer, Integer> singleW = createSingleW(super_node_id);
        for (Integer u : singleW.keySet()) {
            hashCode.code[0] += 1;
            hashCode.code[1] += singleW.get(u) / (superNodeSize * 1.0);
        }
        return hashCode;
    }


    // 生成的HashCode长度太长了，太多0的出现因为度数一般比较小
    public HashCode generateHashCode(int super_node_id) {
        int bin_size = 32; // Integer最大是 (2^32 - 1)
        int k_bins = n / bin_size;
        if (n % bin_size != 0) { k_bins = k_bins + 1; }
        HashCode hashCode = new HashCode(k_bins);
        Map<Integer, Integer> singleW = createSingleW(super_node_id);
        for (Integer u : singleW.keySet()) {
            int bin = u / bin_size;
            hashCode.code[bin] += Math.pow(2, u % bin_size); // 将bin里的第 u%bin_size 位设为1
        }
        return hashCode;
    }

    public Integer generateCode(int super_node_id) {
        Map<Integer, Integer> singleW = createSingleW(super_node_id);
        int prime = 31;
        long code = 1;
        for (Integer u : singleW.keySet()) {
            code = (code * prime +  (long) u * singleW.get(u)) % Integer.MAX_VALUE;
        }
        return (int)code;
    }

    public void originTest2(int iteration, int print_iteration_offset) {
        System.out.println("----------------------------------- ORIGIN ALGORITHM ---------------------------------------");
        int maxGroupSize = 1000;
        Map<Integer, Map<Integer, Record>> it2Records = new HashMap<>();
        for (int it = 1; it <= iteration; it++) {
            List<List<Integer>> groups = new ArrayList<>();
            System.out.println("\n------------------------- ITERATION " + it);
            // DividePhase =============================================================================================
            {
                System.out.println("# Divide Phase");
                long divideStartTime = System.currentTimeMillis();
                int k_bins = signatureLength;
                if (n % k_bins != 0) {
                    k_bins = k_bins + 1;
                }
                // 首先生成长度为k_bins的一个数组用于辅助计算hash签名值
                int[] rot_direction = new int[k_bins];
                Random random = new Random();
                for (int i = 0; i < k_bins; i++) {
                    if (random.nextBoolean()) {
                        rot_direction[i] = 1;
                    } else {
                        rot_direction[i] = -1;
                    }
                }
                // 接着对顶点进行一个编号重排
                randomPermutation();
                // 创建F_OPH数组, 用于存储每个顶点的LSH签名
                OnePermHashSig[] F_OPH_temp = new OnePermHashSig[n];
                for (int A = 0; A < n; A++) {
                    F_OPH_temp[A] = generateSignature(signatureLength, A, rot_direction);
                }
                // 对分组进行排序
                Integer[] G_temp = new Integer[n];
                for (int i = 0; i < n; i++) G_temp[i] = i;
                Arrays.sort(G_temp, (o1, o2) -> OnePermHashSig.compare(F_OPH_temp[o1], F_OPH_temp[o2]));
                // 找到第一个组的开始index
                int g_start_temp = 0;
                while (F_OPH_temp[G_temp[g_start_temp]].unassigned())
                    g_start_temp++;
                // 计算每个组包含哪些超点
                OnePermHashSig g = new OnePermHashSig(k_bins);
                ArrayList<Integer> Q = new ArrayList<>();

                for (int i = g_start_temp; i < n; i++) {
                    // 如果遍历到当前顶点的LSH签名与前面的签名不一样，说明前面组的所有顶点都被识别出来，需要增加新的一组了
                    if (!F_OPH_temp[G_temp[i]].equals(g)) {
                        // 只返回顶点数量大于1的组
                        if (Q.size() > 1){
                            if (Q.size() <= maxGroupSize) {
                                groups.add(Q);
                            } else {
                                List<List<Integer>> subGroups = hierarchicalDivide(Q);
                                System.out.print("Origin group length:" + Q.size() + " divide into group:");
                                for (List<Integer> subGroup : subGroups) {
                                    System.out.print(subGroup.size() + ",");
                                    if(subGroup.size() > 1) groups.add(subGroup);
                                }
                                System.out.println();
                            }
                        }
                        g = F_OPH_temp[G_temp[i]];
                        Q = new ArrayList<>();
                    }
                    Q.add(G_temp[i]);
                }
                if(Q.size() > 1) groups.add(Q);
                System.out.printf("@Time: %5f seconds%n", (System.currentTimeMillis() - divideStartTime) / 1000.0);
            }
            // End =====================================================================================================
            // MergePhase ==============================================================================================
            {
                double threshold = 1 / ((it + 1) * 1.0);
                System.out.println("# Merge Phase");
                System.out.printf("Threshold=%5f%n", threshold);
                // 记录mergePhase的运行时间
                long mergeStartTime = System.currentTimeMillis();
                int total_test = 0;
                int code0_equals = 0;
                int code0_1_equals = 0;
                int code0_1_not = 0;
                int code0_not = 0;
                int merged_num = 0;
                int fails = 0;
                int without_jac = 0;
                // 开始遍历每个组
                Map<Integer, Record> recordMap = new HashMap<>();
                for (List<Integer> group : groups) {
                    long start = System.currentTimeMillis();
                    int success = 0;
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
                        total_test += 1;
                        int counted_jaccard_similarity = 0;
                        for (int j = 0; j < initial_size; j++) {
                            if (hm.get(j) == null)
                                continue;
                            if (j == A) continue;
                            double jaccard_similarity = computeJacSim(hm.get(A), hm.get(j));
                            counted_jaccard_similarity += 1;
                            if (jaccard_similarity > max) {
                                max = jaccard_similarity;
                                idx = j;
                            }
                        }
                        if (idx == -1) {
                            without_jac += 1;
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
                            HashCodeAdjust hashCodeAdjust_1 = generateHashCode_3(group.get(A));
                            HashCodeAdjust hashCodeAdjust_2 = generateHashCode_3(group.get(idx));
                            merged_num += 1;
                            if(hashCodeAdjust_1.code[0]==hashCodeAdjust_2.code[0]){
                                code0_equals += 1;
                                if(hashCodeAdjust_1.code[1] == hashCodeAdjust_2.code[1]){
                                    code0_1_equals += 1;
                                }else{
                                    code0_1_not += 1;
                                }
                            }else{
                                code0_not += 1;
                            }
                            HashMap<Integer, Integer> w_update = updateW(hm.get(A), hm.get(idx));
                            hm.replace(A, w_update);
                            hm.remove(idx);
                            updateSuperNode(group.get(A), group.get(idx));
                            success += 1;
                        } else {
                            fails += 1;
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
                System.out.print("total_test:" + total_test);
                System.out.print(", merged_num:" + merged_num + "(" + (merged_num*1.0/total_test) + ")");
                System.out.print(", code0_equal:" + code0_equals + "(" + (code0_equals*1.0/merged_num) + ")");
                System.out.print(", code0_1_equal:" + code0_1_equals + "(" + (code0_1_equals*1.0/merged_num) + ")");
                System.out.print(", code0_1_not:" + code0_1_not + "(" + (code0_1_not*1.0/merged_num) + ")");
                System.out.print(", code0_not:" + code0_not + "(" + (code0_not*1.0/merged_num) + ")");
                System.out.print(", jaccard similarity fails:" + fails + "(" + (fails*1.0/total_test) + ")");
                System.out.print(", without jaccard similarity:" + without_jac + "(" + (without_jac * 1.0 / total_test) + " )");
                System.out.println();
                System.out.printf("@Time: %5f seconds%n", (System.currentTimeMillis() - mergeStartTime) / 1000.0);
            }
            // End =====================================================================================================
            if (it % print_iteration_offset == 0) {
                System.out.printf("@Time: %5f seconds%n", encodePhase_new());
                evaluatePhase();
            }
        }
    }

    /**
     * 记录origin方法中小组的大小分布和合并成功次数
     * @param iteration 迭代次数
     * @param print_iteration_offset 每隔多少次迭代就进行Encode一次查看压缩率
     */
    public void originTest3(int iteration, int print_iteration_offset) {
        System.out.println("----------------------------------- ORIGIN ALGORITHM ---------------------------------------");
        Map<Integer, Map<Integer, Record>> it2Records = new HashMap<>();
        for (int it = 1; it <= iteration; it++) {
            List<List<Integer>> groups = new ArrayList<>();
            System.out.println("\n------------------------- ITERATION " + it);
            // DividePhase =============================================================================================
            {
                System.out.println("# Divide Phase");
                long divideStartTime = System.currentTimeMillis();
                int k_bins = signatureLength;
                if (n % k_bins != 0) {
                    k_bins = k_bins + 1;
                }
                // 首先生成长度为k_bins的一个数组用于辅助计算hash签名值
                int[] rot_direction = new int[k_bins];
                Random random = new Random();
                for (int i = 0; i < k_bins; i++) {
                    if (random.nextBoolean()) {
                        rot_direction[i] = 1;
                    } else {
                        rot_direction[i] = -1;
                    }
                }
                // 接着对顶点进行一个编号重排
                randomPermutation();
                // 创建F_OPH数组, 用于存储每个顶点的LSH签名
                OnePermHashSig[] F_OPH_temp = new OnePermHashSig[n];
                for (int A = 0; A < n; A++) {
                    F_OPH_temp[A] = generateSignature(signatureLength, A, rot_direction);
                }
                // 对分组进行排序
                Integer[] G_temp = new Integer[n];
                for (int i = 0; i < n; i++) G_temp[i] = i;
                Arrays.sort(G_temp, (o1, o2) -> OnePermHashSig.compare(F_OPH_temp[o1], F_OPH_temp[o2]));
                // 找到第一个组的开始index
                int g_start_temp = 0;
                while (F_OPH_temp[G_temp[g_start_temp]].unassigned())
                    g_start_temp++;
                // 计算每个组包含哪些超点
                OnePermHashSig g = new OnePermHashSig(k_bins);
                ArrayList<Integer> Q = new ArrayList<>();
                for (int i = g_start_temp; i < n; i++) {
                    // 如果遍历到当前顶点的LSH签名与前面的签名不一样，说明前面组的所有顶点都被识别出来，需要增加新的一组了
                    if (!F_OPH_temp[G_temp[i]].equals(g)) {
                        // 只返回顶点数量大于1的组
                        if (Q.size() > 1) groups.add(Q);
                        g = F_OPH_temp[G_temp[i]];
                        Q = new ArrayList<>();
                    }
                    Q.add(G_temp[i]);
                }
                if(Q.size() > 1) groups.add(Q);
                System.out.printf("@Time: %5f seconds%n", (System.currentTimeMillis() - divideStartTime) / 1000.0);
            }
            // End =====================================================================================================

            // MergePhase ==============================================================================================
            {
                double threshold = 1 / ((it + 1) * 1.0);
                System.out.println("# Merge Phase");
                System.out.printf("Threshold=%5f%n", threshold);
                // 记录mergePhase的运行时间
                long mergeStartTime = System.currentTimeMillis();
                int total_test = 0;
                int code0_equals = 0;
                int code0_1_equals = 0;
                int code0_1_not = 0;
                int code0_not = 0;
                int merged_num = 0;
                int fails = 0;
                int without_jac = 0;
                // 开始遍历每个组
                Map<Integer, Record> recordMap = new HashMap<>();
                for (List<Integer> group : groups) {
                    long start = System.currentTimeMillis();
                    int success = 0;
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
                        total_test += 1;
                        int counted_jaccard_similarity = 0;
                        for (int j = 0; j < initial_size; j++) {
                            if (hm.get(j) == null)
                                continue;
                            if (j == A) continue;
                            double jaccard_similarity = computeJacSim(hm.get(A), hm.get(j));
                            counted_jaccard_similarity += 1;
                            if (jaccard_similarity > max) {
                                max = jaccard_similarity;
                                idx = j;
                            }
                        }
                        if (idx == -1) {
                            without_jac += 1;
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
                            HashCodeAdjust hashCodeAdjust_1 = generateHashCode_3(group.get(A));
                            HashCodeAdjust hashCodeAdjust_2 = generateHashCode_3(group.get(idx));
                            merged_num += 1;
                            if(hashCodeAdjust_1.code[0]==hashCodeAdjust_2.code[0]){
                                code0_equals += 1;
                                if(hashCodeAdjust_1.code[1] == hashCodeAdjust_2.code[1]){
                                    code0_1_equals += 1;
                                }else{
                                    code0_1_not += 1;
                                }
                            }else{
                                code0_not += 1;
                            }
//                            if(it > 2 & group.size() > 1000){
//                                System.out.println("===============================================================");
//                                System.out.println("Merge superNode " + group.get(A) + " and " + group.get(idx));
//                                HashCodeAdjust hashCode = generateHashCode_3(group.get(A));
//                                System.out.println(group.get(A) + "'s hasCode: " + hashCode + " neighbors:");
//                                int[] temp = Gr.successorArray(group.get(A));
//                                for (int a : temp) {
//                                    System.out.print(a + ",");
//                                }
//                                System.out.println();
//                                hashCode = generateHashCode_3(group.get(idx));
//                                System.out.println(group.get(idx) + "'s hasCode: " + hashCode + " neighbors:");
//                                temp = Gr.successorArray(group.get(idx));
//                                for (int a : temp) {
//                                    System.out.print(a + ",");
//                                }
//                                System.out.println();
//                                System.out.println("group.size():"+ group.size() + " jaccard_similarity compute times:" + counted_jaccard_similarity);
//                                System.out.println("===============================================================");
//                            }
                            HashMap<Integer, Integer> w_update = updateW(hm.get(A), hm.get(idx));
                            hm.replace(A, w_update);
                            hm.remove(idx);
                            updateSuperNode(group.get(A), group.get(idx));
                            success += 1;
                        } else {
                            fails += 1;
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
                System.out.print("total_test:" + total_test);
                System.out.print(", merged_num:" + merged_num + "(" + (merged_num*1.0/total_test) + ")");
                System.out.print(", code0_equal:" + code0_equals + "(" + (code0_equals*1.0/merged_num) + ")");
                System.out.print(", code0_1_equal:" + code0_1_equals + "(" + (code0_1_equals*1.0/merged_num) + ")");
                System.out.print(", code0_1_not:" + code0_1_not + "(" + (code0_1_not*1.0/merged_num) + ")");
                System.out.print(", code0_not:" + code0_not + "(" + (code0_not*1.0/merged_num) + ")");
                System.out.print(", jaccard similarity fails:" + fails + "(" + (fails*1.0/total_test) + ")");
                System.out.print(", without jaccard similarity:" + without_jac + "(" + (without_jac * 1.0 / total_test) + " )");
                System.out.println();
                System.out.printf("@Time: %5f seconds%n", (System.currentTimeMillis() - mergeStartTime) / 1000.0);
            }
            // End =====================================================================================================
            if (it % print_iteration_offset == 0) {
                System.out.printf("@Time: %5f seconds%n", encodePhase_new());
                evaluatePhase();
            }
        }
        outputResult(it2Records, "origin", 1);
    }



    public void originTest4(int iteration, int print_iteration_offset) {
        System.out.println("----------------------------------- ORIGIN ALGORITHM ---------------------------------------");
        Map<Integer, Map<Integer, Record>> it2Records = new HashMap<>();
        for (int it = 1; it <= iteration; it++) {
            List<List<Integer>> groups = new ArrayList<>();
            System.out.println("\n------------------------- ITERATION " + it);
            // DividePhase =============================================================================================
            {
                System.out.println("# Divide Phase");
                long divideStartTime = System.currentTimeMillis();
                int k_bins = signatureLength;
                if (n % k_bins != 0) {
                    k_bins = k_bins + 1;
                }
                // 首先生成长度为k_bins的一个数组用于辅助计算hash签名值
                int[] rot_direction = new int[k_bins];
                Random random = new Random();
                for (int i = 0; i < k_bins; i++) {
                    if (random.nextBoolean()) {
                        rot_direction[i] = 1;
                    } else {
                        rot_direction[i] = -1;
                    }
                }
                // 接着对顶点进行一个编号重排
                randomPermutation();
                // 创建F_OPH数组, 用于存储每个顶点的LSH签名
                OnePermHashSig[] F_OPH_temp = new OnePermHashSig[n];
                for (int A = 0; A < n; A++) {
                    F_OPH_temp[A] = generateSignature(signatureLength, A, rot_direction);
                }
                // 对分组进行排序
                Integer[] G_temp = new Integer[n];
                for (int i = 0; i < n; i++) G_temp[i] = i;
                Arrays.sort(G_temp, (o1, o2) -> OnePermHashSig.compare(F_OPH_temp[o1], F_OPH_temp[o2]));
                // 找到第一个组的开始index
                int g_start_temp = 0;
                while (F_OPH_temp[G_temp[g_start_temp]].unassigned())
                    g_start_temp++;
                // 计算每个组包含哪些超点
                OnePermHashSig g = new OnePermHashSig(k_bins);
                ArrayList<Integer> Q = new ArrayList<>();
                for (int i = g_start_temp; i < n; i++) {
                    // 如果遍历到当前顶点的LSH签名与前面的签名不一样，说明前面组的所有顶点都被识别出来，需要增加新的一组了
                    if (!F_OPH_temp[G_temp[i]].equals(g)) {
                        // 只返回顶点数量大于1的组
                        if (Q.size() > 1) groups.add(Q);
                        g = F_OPH_temp[G_temp[i]];
                        Q = new ArrayList<>();
                    }
                    Q.add(G_temp[i]);
                }
                if(Q.size() > 1) groups.add(Q);
                System.out.printf("@Time: %5f seconds%n", (System.currentTimeMillis() - divideStartTime) / 1000.0);
            }
            // End =====================================================================================================

            // MergePhase ==============================================================================================
            {
                double threshold = 1 / ((it + 1) * 1.0);
                System.out.println("# Merge Phase");
                System.out.printf("Threshold=%5f%n", threshold);
                // 记录mergePhase的运行时间
                long mergeStartTime = System.currentTimeMillis();
                // 开始遍历每个组
                Map<Integer, Record> recordMap = new HashMap<>();
                for (List<Integer> group : groups) {
                    long start = System.currentTimeMillis();
                    int success = 0;
                    HashMap<Integer, HashMap<Integer, Integer>> hm = createW(group);
                    int initial_size = hm.size();
                    boolean flag = group.size() > 1000 && it >= 3;
                    while (hm.size() > 1) {
                        Random rand = new Random();
                        // 从组内随机找到一个超点A
                        int A = rand.nextInt(initial_size);
                        if (hm.get(A) == null) continue;
                        // 变量max记录最大的Jaccard_similarity
                        double max = 0;
                        double max_saving = 0;
                        double savings = 0;
                        // 变量idx记录最大Jaccard_similarity的那个顶点
                        int idx = -1;
                        int idx_saving = -1;
                        // 遍历组内其他顶点，找到与A的Jaccard Similarity最大的那个顶点
                        for (int j = 0; j < initial_size; j++) {
                            if (hm.get(j) == null)
                                continue;
                            if (j == A) continue;
                            double jaccard_similarity = computeJacSim(hm.get(A), hm.get(j));
                            if(flag){
                                savings = computeSaving(hm.get(A), hm.get(j), group.get(A), group.get(j));
                                System.out.print(group.get(j) + ", " + jaccard_similarity + ", " + savings);
                                Map<Integer, Integer> node2weights = new HashMap<>();
                                for(Integer u : hm.get(j).keySet()){
                                    int edges = node2weights.getOrDefault(S[u], 0);
                                    node2weights.put(S[u], edges + hm.get(j).get(u));
                                }
                                for (Integer u : node2weights.keySet()) {
                                    System.out.print(", (" + u + "," + (node2weights.get(u)*1.0) + ")");
                                }
                                System.out.println();
                            }
                            if (jaccard_similarity > max) {
                                max = jaccard_similarity;
                                idx = j;
                            }
                            if(flag && savings > max_saving){
                                max_saving = savings;
                                idx_saving = j;
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
                        savings = computeSaving(hm.get(A), hm.get(idx), group.get(A), group.get(idx));
                        if(flag){
                            System.out.println("max_jaccard_index:" + idx + ":" + max + ":" + savings +", max_saving_index:" + idx_saving + ":" + max_saving  + "======");
                            System.out.println();

                            System.out.println("== superNode A ==========================================");
                            int[] nodes = recoverSuperNode(group.get(A));
                            System.out.print("id: " + group.get(A) + ", node_size: " + nodes.length + "(");
                            for(int u : nodes){
                                System.out.print(u + ",");
                            }
                            System.out.println(")");
                            System.out.print("neighbors contain: ");
                            Map<Integer, Integer> candidate_size = new HashMap<>();
                            Map<Integer, Integer> edges_num = new HashMap<>();
                            for(Integer u : hm.get(A).keySet()){
                                System.out.print("(" + u + "," + hm.get(A).get(u) + "), ");
                                if(!candidate_size.containsKey(S[u])){
                                    int[] Nodes = recoverSuperNode(S[u]);
                                    candidate_size.put(S[u], Nodes.length);
                                    edges_num.put(S[u], hm.get(A).get(u));
                                }else{
                                    edges_num.put(S[u], edges_num.get(S[u]) + hm.get(A).get(u));
                                }
                            }
                            System.out.print("\n");
                            double cost = 0;
                            for (Integer u : edges_num.keySet()) {
                                int E = edges_num.get(u);
                                double compare = u == group.get(A) ? ((nodes.length * 1.0 * (nodes.length - 1)) / 2.0) : (nodes.length * 1.0 * candidate_size.get(u));
                                cost = (E <= compare / 2.0) ? (E) : (1 + compare - E);
                                System.out.print("id: " + u + ", cost: " + cost +", actual_edges: " + E + ", node_size: " + candidate_size.get(u) + "(");
                                nodes = recoverSuperNode(u);
                                for(int v : nodes){
                                    System.out.print(v + ",");
                                }
                                System.out.println(")");
                            }
                            System.out.println("=========================================================\n");

                            System.out.println("== superNode B with max jaccard =============================");
                            nodes = recoverSuperNode(group.get(idx));
                            System.out.print("id: " + group.get(idx) + ", node_size: " + nodes.length + "(");
                            for(int u : nodes){
                                System.out.print(u + ",");
                            }
                            System.out.println(")");
                            System.out.println("jaccard:" + max + ", saving:" + savings);
                            System.out.print("neighbors contain: ");
                            candidate_size = new HashMap<>();
                            edges_num = new HashMap<>();
                            for(Integer u : hm.get(idx).keySet()){
                                System.out.print("(" + u + "," + hm.get(idx).get(u) + "), ");
                                if(!candidate_size.containsKey(S[u])){
                                    int[] Nodes = recoverSuperNode(S[u]);
                                    candidate_size.put(S[u], Nodes.length);
                                    edges_num.put(S[u], hm.get(idx).get(u));
                                }else{
                                    edges_num.put(S[u], edges_num.get(S[u]) + hm.get(idx).get(u));
                                }
                            }
                            System.out.println();
                            for (Integer u : edges_num.keySet()) {
                                int E = edges_num.get(u);
                                double compare = u == group.get(idx) ? ((nodes.length * 1.0 * (nodes.length - 1)) / 2.0) : (nodes.length * 1.0 * candidate_size.get(u));
                                cost = (E <= compare / 2.0) ? (E) : (1 + compare - E);
                                System.out.print("id: " + u + ", cost: " + cost +", actual_edges: " + E + ", node_size: " + candidate_size.get(u) + "(");
                                nodes = recoverSuperNode(u);
                                for(int v : nodes){
                                    System.out.print(v + ",");
                                }
                                System.out.println(")");
                            }
                            System.out.println("=========================================================\n");

                            System.out.println("== superNode B with max saving =============================");
                            double max_jaccard = computeJacSim(hm.get(A), hm.get(idx_saving));
                            nodes = recoverSuperNode(group.get(idx_saving));
                            System.out.print("id: " + group.get(idx_saving) + ", node_size: " + nodes.length + "(");
                            for(int u : nodes){
                                System.out.print(u + ",");
                            }
                            System.out.println(")");
                            System.out.println("jaccard:" + max_jaccard + ", saving:" + max_saving);
                            System.out.print("neighbors contain: ");
                            candidate_size = new HashMap<>();
                            edges_num = new HashMap<>();
                            for(Integer u : hm.get(idx_saving).keySet()){
                                System.out.print("(" + u + "," + hm.get(idx_saving).get(u) + "), ");
                                if(!candidate_size.containsKey(S[u])){
                                    int[] Nodes = recoverSuperNode(S[u]);
                                    candidate_size.put(S[u], Nodes.length);
                                    edges_num.put(S[u], hm.get(idx_saving).get(u));
                                }else{
                                    edges_num.put(S[u], edges_num.get(S[u]) + hm.get(idx_saving).get(u));
                                }
                            }
                            System.out.println();
                            for (Integer u : edges_num.keySet()) {
                                int E = edges_num.get(u);
                                double compare = u == group.get(idx_saving) ? ((nodes.length * 1.0 * (nodes.length - 1)) / 2.0) : (nodes.length * 1.0 * candidate_size.get(u));
                                cost = (E <= compare / 2.0) ? (E) : (1 + compare - E);
                                System.out.print("id: " + u + ", cost: " + cost +", actual_edges: " + E + ", node_size: " + candidate_size.get(u) + "(");
                                nodes = recoverSuperNode(u);
                                for(int v : nodes){
                                    System.out.print(v + ",");
                                }
                                System.out.println(")");
                            }
                            System.out.println("=========================================================\n");
                            break;
                        }
                        if (savings >= threshold) {
//                            if(it >= 5){
//                                System.out.println("Merge " + group.get(A) + " and " + group.get(idx) + ", jac=" + max + ", saving=" + savings);
//                                int A_common = 0;
//                                int A_different = 0;
//                                int B_common = 0;
//                                int B_different = 0;
//                                int nums = 0;
//                                Map<Integer, Integer> common_neighbors = new HashMap<>();
//                                int A_length = superNodeLength(group.get(A));
//                                int B_length = superNodeLength(group.get(idx));
//                                Map<Integer, Integer> A_own = new HashMap<>();
//                                Map<Integer, Integer> B_own = new HashMap<>();
//                                for(Integer u : hm.get(A).keySet()){
//                                    int origin = 0;
//                                    if(A_own.containsKey(S[u])){
//                                        origin = A_own.get(S[u]);
//                                    }
//                                    A_own.put(S[u], origin + hm.get(A).get(u));
//                                }
//                                for (Integer u : hm.get(idx).keySet()) {
//                                    int origin = 0;
//                                    if (B_own.containsKey(S[u])) {
//                                        origin = B_own.get(S[u]);
//                                    }
//                                    B_own.put(S[u], origin + hm.get(idx).get(u));
//                                }
//                                for(Integer u : A_own.keySet()){
//                                    if (B_own.containsKey(u)) {
//                                        int length = superNodeLength(u);
//                                        System.out.print("(" + ((A_own.get(u)*1.0)/(length*A_length)) + "," + ((B_own.get(u)*1.0)/(length*B_length)) + "), ");
//                                    }
//                                }
//                                System.out.println();
////                                System.out.println("common:" + nums + "(A:" + A_common + ",B:"+B_common+")"+", different_A:" + A_different + ", different_B:"+B_different);
//                            }
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
                System.out.printf("@Time: %5f seconds%n", (System.currentTimeMillis() - mergeStartTime) / 1000.0);
            }
            // End =====================================================================================================
            if (it % print_iteration_offset == 0) {
                System.out.printf("@Time: %5f seconds%n", encodePhase_new());
                evaluatePhase();
            }
        }
    }


    public void originTest_outputResults(int iteration, int print_iteration_offset) {
        System.out.println("----------------------------------- ORIGIN ALGORITHM ---------------------------------------");
        Map<Integer, Map<Integer, Record>> it2Records = new HashMap<>();
        for (int it = 1; it <= iteration; it++) {
            List<List<Integer>> groups = new ArrayList<>();
            System.out.println("\n------------------------- ITERATION " + it);
            // DividePhase =============================================================================================
            {
                System.out.println("# Divide Phase");
                long divideStartTime = System.currentTimeMillis();
                int k_bins = signatureLength;
                if (n % k_bins != 0) {
                    k_bins = k_bins + 1;
                }
                // 首先生成长度为k_bins的一个数组用于辅助计算hash签名值
                int[] rot_direction = new int[k_bins];
                Random random = new Random();
                for (int i = 0; i < k_bins; i++) {
                    if (random.nextBoolean()) {
                        rot_direction[i] = 1;
                    } else {
                        rot_direction[i] = -1;
                    }
                }
                // 接着对顶点进行一个编号重排
                randomPermutation();
                // 创建F_OPH数组, 用于存储每个顶点的LSH签名
                OnePermHashSig[] F_OPH_temp = new OnePermHashSig[n];
                for (int A = 0; A < n; A++) {
                    F_OPH_temp[A] = generateSignature(signatureLength, A, rot_direction);
                }
                // 对分组进行排序
                Integer[] G_temp = new Integer[n];
                for (int i = 0; i < n; i++) G_temp[i] = i;
                Arrays.sort(G_temp, (o1, o2) -> OnePermHashSig.compare(F_OPH_temp[o1], F_OPH_temp[o2]));
                // 找到第一个组的开始index
                int g_start_temp = 0;
                while (F_OPH_temp[G_temp[g_start_temp]].unassigned())
                    g_start_temp++;
                // 计算每个组包含哪些超点
                OnePermHashSig g = new OnePermHashSig(k_bins);
                ArrayList<Integer> Q = new ArrayList<>();
                for (int i = g_start_temp; i < n; i++) {
                    // 如果遍历到当前顶点的LSH签名与前面的签名不一样，说明前面组的所有顶点都被识别出来，需要增加新的一组了
                    if (!F_OPH_temp[G_temp[i]].equals(g)) {
                        // 只返回顶点数量大于1的组
                        if (Q.size() > 1) groups.add(Q);
                        g = F_OPH_temp[G_temp[i]];
                        Q = new ArrayList<>();
                    }
                    Q.add(G_temp[i]);
                }
                if(Q.size() > 1) groups.add(Q);
                System.out.printf("@Time: %5f seconds%n", (System.currentTimeMillis() - divideStartTime) / 1000.0);
            }
            // End =====================================================================================================

            // MergePhase ==============================================================================================
            {
                double threshold = 1 / ((it + 1) * 1.0);
                System.out.println("# Merge Phase");
                System.out.printf("Threshold=%5f%n", threshold);
                // 记录mergePhase的运行时间
                long mergeStartTime = System.currentTimeMillis();
                // 开始遍历每个组
                Map<Integer, Record> recordMap = new HashMap<>();
                for (List<Integer> group : groups) {
                    long start = System.currentTimeMillis();
                    int success = 0;
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
                System.out.printf("@Time: %5f seconds%n", (System.currentTimeMillis() - mergeStartTime) / 1000.0);
            }
            // End =====================================================================================================

            String filename = "/home/hallo/下载/GraphSummarization/result/20211212/superNode_" + it + ".csv";
            File writeFile = new File(filename);
            try{
                BufferedWriter out = new BufferedWriter(new FileWriter(writeFile));
                out.write("node,S,I,J\n");
                out.flush();
                for (int i = 0; i < n; i++) {
                    out.write(i + "," + S[i] + ","  + I[i] + "," + J[i] +"\n");
                    out.flush();
                }
                out.close();
            }catch (Exception e){
                System.out.println(e.getMessage());
            }
            System.out.println("Create file " + filename + " success!");
            if (it % print_iteration_offset == 0) {
                System.out.printf("@Time: %5f seconds%n", encodePhase_new());
                evaluatePhase();
            }
        }
//        outputResult(it2Records, "origin", 1);
    }


    public void originTest_recoverSuperNode(int iteration){
        // Recover results =============================================================================================
        String filename = "/home/hallo/下载/GraphSummarization/result/20211212/superNode_" + iteration + ".csv";
        File readFile = new File(filename);
        try{
            BufferedReader in = new BufferedReader(new FileReader(readFile));
            in.readLine();
            String line;
            while ((line = in.readLine()) != null) {
                int node = Integer.parseInt(line.split(",")[0]);
                int S_node = Integer.parseInt(line.split(",")[1]);
                int I_node = Integer.parseInt(line.split(",")[2]);
                int J_node = Integer.parseInt(line.split(",")[3]);
                S[node] = S_node;
                I[node] = I_node;
                J[node] = J_node;
            }
            in.close();
        }catch (Exception e){
            System.out.println(e.getMessage());
        }
        // End =========================================================================================================

        List<List<Integer>> groups = new ArrayList<>();
        OnePermHashSig[] F_OPH_temp = new OnePermHashSig[n];
        // Divide Phase ================================================================================================
        {
            System.out.println("# Divide Phase");
            long divideStartTime = System.currentTimeMillis();
            int k_bins = signatureLength;
            if (n % k_bins != 0) {
                k_bins = k_bins + 1;
            }
            // 首先生成长度为k_bins的一个数组用于辅助计算hash签名值
            int[] rot_direction = new int[k_bins];
            Random random = new Random();
            for (int i = 0; i < k_bins; i++) {
                if (random.nextBoolean()) {
                    rot_direction[i] = 1;
                } else {
                    rot_direction[i] = -1;
                }
            }
            // 接着对顶点进行一个编号重排
            randomPermutation();
            // 创建F_OPH数组, 用于存储每个顶点的LSH签名
//            OnePermHashSig[] F_OPH_temp = new OnePermHashSig[n];
            for (int A = 0; A < n; A++) {
                F_OPH_temp[A] = generateSignature(signatureLength, A, rot_direction);
            }
            // 对分组进行排序
            Integer[] G_temp = new Integer[n];
            for (int i = 0; i < n; i++) G_temp[i] = i;
            Arrays.sort(G_temp, (o1, o2) -> OnePermHashSig.compare(F_OPH_temp[o1], F_OPH_temp[o2]));
            // 找到第一个组的开始index
            int g_start_temp = 0;
            while (F_OPH_temp[G_temp[g_start_temp]].unassigned())
                g_start_temp++;
            // 计算每个组包含哪些超点
            OnePermHashSig g = new OnePermHashSig(k_bins);
            ArrayList<Integer> Q = new ArrayList<>();
            for (int i = g_start_temp; i < n; i++) {
                // 如果遍历到当前顶点的LSH签名与前面的签名不一样，说明前面组的所有顶点都被识别出来，需要增加新的一组了
                if (!F_OPH_temp[G_temp[i]].equals(g)) {
                    // 只返回顶点数量大于1的组
                    if (Q.size() > 1) groups.add(Q);
                    g = F_OPH_temp[G_temp[i]];
                    Q = new ArrayList<>();
                }
                Q.add(G_temp[i]);
            }
            if(Q.size() > 1) groups.add(Q);
            System.out.printf("@Time: %5f seconds%n", (System.currentTimeMillis() - divideStartTime) / 1000.0);
        }
        // End =========================================================================================================

        int[] subGroup = new int[]{108738, 108741, 108742, 108751, 108762, 108768, 108775, 108778, 110462, 110465, 110549};
        Map<Integer, HashMap<Integer,Integer>> hm = createW(subGroup, subGroup.length);
        Set<Integer> min_common_neighbors = new TreeSet<>();
        Set<Integer> max_common_neighbors = new TreeSet<>();
        for(int i=0; i<subGroup.length; i++){
            System.out.println("node:" + subGroup[i] + ", signature:" + F_OPH_temp[subGroup[i]].toString());
            Map<Integer, Integer> candidate_superNode = new HashMap<>();
            Map<Integer, Integer> candidate_superNode_edges = new HashMap<>();
            System.out.print("neighbors contain:");
            for (Integer u : hm.get(i).keySet()) {
                int[] candidate_superNode_backup = recoverSuperNode(S[u]);
                if(!candidate_superNode.containsKey(S[u])){
                    candidate_superNode.put(S[u], candidate_superNode_backup.length);
                    candidate_superNode_edges.put(S[u], hm.get(i).get(u));
                }else{
                    candidate_superNode_edges.put(S[u], candidate_superNode_edges.get(S[u]) + hm.get(i).get(u));
                }
            }
            for (Integer u : candidate_superNode_edges.keySet()) {
                int[] A = recoverSuperNode(S[subGroup[i]]);
                int E = candidate_superNode_edges.get(u);
                double compare = u == subGroup[i] ? ((A.length*1.0*(A.length-1))/2.0) : (A.length*1.0*candidate_superNode.get(u));
                System.out.print("(" + A.length + "," + u + "," + (E/compare) + "), ");
            }
            System.out.println();
            for(int j=i+1; j<subGroup.length; j++){
                double jaccard = computeJacSim(hm.get(i), hm.get(j));
                double saving = computeSaving(hm.get(i), hm.get(j), subGroup[i], subGroup[j]);
                System.out.print("(" + subGroup[j] + "," + jaccard + "," + saving + ") ");
            }
            System.out.println("\n=============================================================================");
        }
//        List<Integer> store_group = new ArrayList<>();
//        List<List<Integer>> store_groups = new ArrayList<>();
//        boolean flag = true;
//        for (List<Integer> group : groups) {
//            if(group.size() > 1000){
//                System.out.print("Group contain " + group.size() + " nodes, divide again:");
//                int max = 0;
//                List<List<Integer>> subGroups = hierarchicalDivide(signatureLength, group);
//                if(flag){
//                    store_group.addAll(group);
//                    store_groups.addAll(subGroups);
//                    flag = false;
//                }
//                for (List<Integer> subGroup : subGroups) {
//                    if(subGroup.size() > max) max = subGroup.size();
//                    System.out.print(subGroup.size() + ",");
//                }
//                System.out.println("\nmax group size:" + max + ", decrease " + ((group.size()-max)/(group.size()*1.0)));
//            }
//        }
//
//        Map<Integer, Integer> node2SubGroup = new HashMap<>();
//        for (int i=0; i< store_groups.size(); i++) {
//            List<Integer> subGroup = store_groups.get(i);
//            for (Integer u : subGroup) {
//                node2SubGroup.put(u, i);
//            }
//        }
//        for (Integer u : store_group) {
//            if(!node2SubGroup.containsKey(u)){
//                node2SubGroup.put(u, -1);
//            }
//        }
//        filename = "/home/hallo/下载/GraphSummarization/result/20211213/groups_" + iteration + ".csv";
//        File outFile = new File(filename);
//        try{
//            BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
//            out.write("node,new_group\n");
//            out.flush();
//            for (Integer u : node2SubGroup.keySet()) {
//                out.write(u + "," + node2SubGroup.get(u)+"\n");
//                out.flush();
//            }
//            out.close();
//        }catch (Exception e){
//            System.out.println(e.getMessage());
//        }
        // Evaluate ====================================================================================================
//        System.out.printf("@Time: %5f seconds%n", encodePhase_new());
//        evaluatePhase();
    }

    /**
     * 记录origin方法中小组的大小分布和合并成功次数
     * @param iteration 迭代次数
     * @param print_iteration_offset 每隔多少次迭代就进行Encode一次查看压缩率
     */
    public void sequentialTest(int iteration, int print_iteration_offset) {
        System.out.println("---------------------------------- SEQUENTIAL ALGORITHM ------------------------------------");
        Map<Integer, Map<Integer, Record>> it2Records = new HashMap<>();
        for (int it = 1; it <= iteration; it++) {
            List<List<Integer>> groups = new ArrayList<>();
            System.out.println("\n------------------------- ITERATION " + it);
            // DividePhase =============================================================================================
            {
                int maxGroupSize = 1000;
                System.out.println("# Divide Phase");
                long divideStartTime = System.currentTimeMillis();
                int k_bins = signatureLength;
                if (n % k_bins != 0) {
                    k_bins = k_bins + 1;
                }
                // 首先生成长度为k_bins的一个数组用于辅助计算hash签名值
                int[] rot_direction = new int[k_bins];
                Random random = new Random();
                for (int i = 0; i < k_bins; i++) {
                    if (random.nextBoolean()) {
                        rot_direction[i] = 1;
                    } else {
                        rot_direction[i] = -1;
                    }
                }
                // 接着对顶点进行一个编号重排
                randomPermutation();
                // 创建F_OPH数组, 用于存储每个顶点的LSH签名
                OnePermHashSig[] F_OPH_temp = new OnePermHashSig[n];
                for (int A = 0; A < n; A++) {
                    F_OPH_temp[A] = generateSignature(signatureLength, A, rot_direction);
                }
                // 对分组进行排序
                Integer[] G_temp = new Integer[n];
                for (int i = 0; i < n; i++) G_temp[i] = i;
                Arrays.sort(G_temp, (o1, o2) -> OnePermHashSig.compare(F_OPH_temp[o1], F_OPH_temp[o2]));
                // 找到第一个组的开始index
                int g_start_temp = 0;
                while (F_OPH_temp[G_temp[g_start_temp]].unassigned())
                    g_start_temp++;
                // 计算每个组包含哪些超点
                OnePermHashSig g = new OnePermHashSig(k_bins);
                ArrayList<Integer> Q = new ArrayList<>();
                for (int i = g_start_temp; i < n; i++) {
                    // 如果遍历到当前顶点的LSH签名与前面的签名不一样，说明前面组的所有顶点都被识别出来，需要增加新的一组了
                    if (!F_OPH_temp[G_temp[i]].equals(g)) {
                        // 只返回顶点数量大于1的组
                        if (Q.size() > 1){
                            if (Q.size() <= maxGroupSize) {
                                groups.add(Q);
                            } else {
                                List<List<Integer>> subGroups = sequentialSplitGroup(Q, maxGroupSize);
                                for (List<Integer> subGroup : subGroups) {
                                    if(subGroup.size() > 1) groups.add(subGroup);
                                }
                            }
                        }
                        g = F_OPH_temp[G_temp[i]];
                        Q = new ArrayList<>();
                    }
                    Q.add(G_temp[i]);
                }
                if(Q.size() > 1) groups.add(Q);
                System.out.printf("@Time: %5f seconds%n", (System.currentTimeMillis() - divideStartTime) / 1000.0);
            }
            // End =====================================================================================================

            // MergePhase ==============================================================================================
            {
                double threshold = 1 / ((it + 1) * 1.0);
                System.out.println("# Merge Phase");
                System.out.printf("Threshold=%5f%n", threshold);
                // 记录mergePhase的运行时间
                long mergeStartTime = System.currentTimeMillis();
                // 开始遍历每个组
                Map<Integer, Record> recordMap = new HashMap<>();
                for (List<Integer> group : groups) {
                    long start = System.currentTimeMillis();
                    int success = 0;
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
                System.out.printf("@Time: %5f seconds%n", (System.currentTimeMillis() - mergeStartTime) / 1000.0);
            }
            // End =====================================================================================================
            if (it % print_iteration_offset == 0) {
                System.out.printf("@Time: %5f seconds%n", encodePhase_new());
                evaluatePhase();
            }
        }
        outputResult(it2Records, "sequential", 1);
    }

    /**
     * 记录origin方法中小组的大小分布和合并成功次数
     * @param iteration 迭代次数
     * @param print_iteration_offset 每隔多少次迭代就进行Encode一次查看压缩率
     */
    public void hierarchicalTest(int iteration, int print_iteration_offset) {
        System.out.println("---------------------------------- HIERARCHICAL ALGORITHM ----------------------------------");
        Map<Integer, Map<Integer, Record>> it2Records = new HashMap<>();
        for (int it = 1; it <= iteration; it++) {
            List<List<Integer>> groups = new ArrayList<>();
            System.out.println("\n------------------------- ITERATION " + it);
            // DividePhase =============================================================================================
            {
                int maxGroupSize = 1000;
                System.out.println("# Divide Phase");
                long divideStartTime = System.currentTimeMillis();
                int k_bins = signatureLength;
                if (n % k_bins != 0) {
                    k_bins = k_bins + 1;
                }
                // 首先生成长度为k_bins的一个数组用于辅助计算hash签名值
                int[] rot_direction = new int[k_bins];
                Random random = new Random();
                for (int i = 0; i < k_bins; i++) {
                    if (random.nextBoolean()) {
                        rot_direction[i] = 1;
                    } else {
                        rot_direction[i] = -1;
                    }
                }
                // 接着对顶点进行一个编号重排
                randomPermutation();
                // 创建F_OPH数组, 用于存储每个顶点的LSH签名
                OnePermHashSig[] F_OPH_temp = new OnePermHashSig[n];
                for (int A = 0; A < n; A++) {
                    F_OPH_temp[A] = generateSignature(signatureLength, A, rot_direction);
                }
                // 对分组进行排序
                Integer[] G_temp = new Integer[n];
                for (int i = 0; i < n; i++) G_temp[i] = i;
                Arrays.sort(G_temp, (o1, o2) -> OnePermHashSig.compare(F_OPH_temp[o1], F_OPH_temp[o2]));
                // 找到第一个组的开始index
                int g_start_temp = 0;
                while (F_OPH_temp[G_temp[g_start_temp]].unassigned())
                    g_start_temp++;
                // 计算每个组包含哪些超点
                OnePermHashSig g = new OnePermHashSig(k_bins);
                ArrayList<Integer> Q = new ArrayList<>();
                int before_max_group = 0;
                int after_max_group = 0;
                for (int i = g_start_temp; i < n; i++) {
                    // 如果遍历到当前顶点的LSH签名与前面的签名不一样，说明前面组的所有顶点都被识别出来，需要增加新的一组了
                    if (!F_OPH_temp[G_temp[i]].equals(g)) {
                        // 只返回顶点数量大于1的组
                        if (Q.size() > 1){
                            if (Q.size() <= maxGroupSize) {
                                groups.add(Q);
                            } else {
                                if(Q.size() > before_max_group) before_max_group = Q.size();
                                System.out.print("group contain " + Q.size() + " nodes divide into more groups:");
                                List<List<Integer>> subGroups = hierarchicalDivide(2*signatureLength, Q);
                                for (List<Integer> subGroup : subGroups) {
                                    System.out.print(subGroup.size() + ",");
                                    if(subGroup.size() > 1) groups.add(subGroup);
                                    if(subGroup.size() > after_max_group) after_max_group = subGroup.size();
                                }
                                System.out.println();
                            }
                        }
                        g = F_OPH_temp[G_temp[i]];
                        Q = new ArrayList<>();
                    }
                    Q.add(G_temp[i]);
                }
                if(Q.size() > 1) {
                    groups.add(Q);
                    if(Q.size() > before_max_group) before_max_group = Q.size();
                }
                System.out.println("Before max group size:" + before_max_group + ", after max group size:" + after_max_group);
                System.out.printf("@Time: %5f seconds%n", (System.currentTimeMillis() - divideStartTime) / 1000.0);
            }
            // End =====================================================================================================

            // MergePhase ==============================================================================================
            {
                double threshold = 1 / ((it + 1) * 1.0);
                System.out.println("# Merge Phase");
                System.out.printf("Threshold=%5f%n", threshold);
                // 记录mergePhase的运行时间
                long mergeStartTime = System.currentTimeMillis();
                // 开始遍历每个组
                Map<Integer, Record> recordMap = new HashMap<>();
                for (List<Integer> group : groups) {
                    long start = System.currentTimeMillis();
                    int success = 0;
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
                System.out.printf("@Time: %5f seconds%n", (System.currentTimeMillis() - mergeStartTime) / 1000.0);
            }
            // End =====================================================================================================
            if (it % print_iteration_offset == 0) {
                System.out.printf("@Time: %5f seconds%n", encodePhase_new());
                evaluatePhase();
            }
        }
        outputResult(it2Records, "hierarchical.txt", 1);
    }

    /**
     * @param iteration              迭代次数
     * @param print_iteration_offset 每执行多少次迭代就进行一次 encode 和 evaluate 进行结果输出
     */
    @Override
    public void run(int iteration, int print_iteration_offset) {
//        originTest_outputResults(iteration, print_iteration_offset);
//        originTest_recoverSuperNode(2);
//        originTest4(iteration, print_iteration_offset);
//        originTest2(iteration, print_iteration_offset);
//        originTest(iteration, print_iteration_offset);
//        sequentialTest(iteration, print_iteration_offset);
//        hierarchicalTest(iteration, print_iteration_offset);
//        System.out.println("----------------------------------- Test ALGORITHM ----------------------------------------");
//        for (int it = 1; it <= iteration; it++) {
//            System.out.println("\n------------------------- ITERATION " + it);
//            double threshold = 1 / ((it + 1) * 1.0);
//            List<List<Integer>> groups = dividePhase_new(signatureLength);
//            System.out.printf("@Time: %5f seconds%n", mergePhase_new(threshold, groups));
//            if (it % print_iteration_offset == 0) {
//                System.out.printf("@Time: %5f seconds%n", encodePhase_new());
//                evaluatePhase();
//            }
//        }
//        System.out.println("----------------------------------- LDME ALGORITHM ----------------------------------------");
        long starTime = System.currentTimeMillis();
        for (int it = 1; it <= iteration; it++) {
            logger_.info(String.format("迭代轮数: %d", it));
            double threshold = 1 / ((it + 1) * 1.0);
            logger_.info(String.format("分组结束, 花费时间 %5f seconds", dividePhase()));
            logger_.info(String.format("合并结束, 花费时间 %5f seconds", mergePhase_test2(it, threshold)));
            if (it % print_iteration_offset == 0) {
                logger_.info(String.format("编码结束, 花费时间 %5f seconds", encodePhase_test()));
                evaluatePhase();
                testRecoverNeighbors(n, "test");
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
        logger_.info(String.format("程序运行结束, 总花费时间 %5f seconds", (System.currentTimeMillis() - starTime) / 1000.0));

//        int i = 1;
//        for (Double m : merged_success) {
//            System.out.println("Iteration " + i + ", average merged success pairs:" + m);
//            i++;
//        }

//        outputResult("/home/hallo/下载/GraphSummarization/result/20211102/LDME_5.txt");
    }
}
