package graph_summarization;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
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

    /**
     * 构造函数，继承了Summary父类的构造函数，同时初始化自己的数据结构
     *
     * @param basename 数据集的基本名字
     * @throws Exception
     */
    public SWeG(String basename) throws Exception {
        super(basename);
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
        int[] neighbors = Gr.successorArray(u);
        for (int i = 0; i < neighbors.length; i++) {
            int v = neighbors[i];
            if (f_u > h[v]) {
                f_u = h[v];
            }
        }
        return f_u;
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
        System.out.println("# Divide Phase");
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
     * @param iteration              迭代次数
     * @param print_iteration_offset 每执行多少次迭代就进行一次 encode 和 evaluate 进行结果输出
     */
    @Override
    public void run(int iteration, int print_iteration_offset) {
        System.out.println("----------------------------------- SWeG ALGORITHM ----------------------------------------");
        for (int it = 1; it <= iteration; it++) {
            System.out.println("\n------------------------- ITERATION " + it);
            double threshold = 1 / ((it + 1) * 1.0);
//            double Threshold = 0.5 - it * 0.05;
            System.out.println(String.format("@Time: %5f seconds", dividePhase()));
            System.out.println(String.format("@Time: %5f seconds", mergePhase(threshold)));
            if (it % print_iteration_offset == 0) {
                System.out.println(String.format("@Time: %5f seconds", encodePhase()));
                evaluatePhase();
            }
        }
    }
}
