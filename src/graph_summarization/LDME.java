package graph_summarization;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

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

    /**
     * 构造函数，用于初始化一些共同的结构
     *
     * @param basename 数据集的基本名字
     * @throws Exception
     */
    public LDME(String basename, int signatureLength) throws Exception {
        super(basename);
        this.signatureLength = signatureLength;
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

        // 初始化F_OPH数组, 用于存储每个顶点的哈希值
        F_OPH = new OnePermHashSig[n];
        for(int A=0; A<n ; A++)
            F_OPH[A] = new OnePermHashSig(k_bins);

        for (int A = 0; A < n; A++) {
            // A不是一个超点
            if (I[A] == -1) continue;
            for (int v = I[A]; ; v=J[v]) {
                int[] neighbours = Gr.successorArray(v);
                for (int j = 0; j < neighbours.length; j++) {
                    int permuted_h = h[neighbours[j]];
                    int permuted_bin = permuted_h / bin_size;
                    if (F_OPH[A].sig[permuted_bin] == -1 || permuted_h % bin_size < F_OPH[A].sig[permuted_bin]) {
                        F_OPH[A].sig[permuted_bin] = permuted_h % bin_size;
                    }
                }
                if(J[v]==-1)
                    break;
            }

            // rotation
            for (int A_bin = 0; A_bin < k_bins; A_bin++) {
                int direction = rot_direction[A_bin];
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

        return (System.currentTimeMillis() - startTime) / 1000.0;
    }

    /**
     * 合并阶段，LDME算法在每个小组内的顶点采用Random方式进行合并
     */
    @Override
    public double mergePhase(double threshold){
        System.out.println("# Merge Phase");
        System.out.println(String.format("Threshold=%5f", threshold));
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
     * @param iteration              迭代次数
     * @param print_iteration_offset 每执行多少次迭代就进行一次 encode 和 evaluate 进行结果输出
     */
    @Override
    public void run(int iteration, int print_iteration_offset) {
        System.out.println("----------------------------------- LDME ALGORITHM ----------------------------------------");
        for (int it = 1; it <= iteration; it++) {
            System.out.println("\n------------------------- ITERATION " + it);
            double threshold = 1 / ((it + 1) * 1.0);
//            double Threshold = 0.5 - it * 0.05;
            System.out.println(String.format("@Time: %5f seconds", dividePhase()));
            System.out.println(String.format("@Time: %5f seconds", mergePhase(threshold)));
            if (it % print_iteration_offset == 0) {
//                System.out.println(String.format("@Time: %5f seconds", encodePhase()));
                System.out.println(String.format("@Time: %5f seconds", encodePhase_new()));
                evaluatePhase();
            }
        }
    }
}
