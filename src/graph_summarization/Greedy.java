package graph_summarization;

import java.util.*;

public class Greedy extends Summary{

    // 记录所有顶点的two-hops邻居集合
    Set<Integer>[] two_hops_neighbors;
    // 优先队列
    Queue<NodesPair> H;
    // 用于得到位于优先队列H中包含某个超点的顶点对
    HashMap<Integer, HashSet<Integer>> H_record;
    // 记录所有顶点的邻居边信息
    HashMap<Integer, HashMap<Integer, Integer>> all_W;

    /**
     * 构造函数，用于初始化一些共同的结构
     *
     * @param basename 数据集的基本名字
     * @throws Exception
     */
    public Greedy(String basename) throws Exception {
        super(basename);
        H = new PriorityQueue<>();
        H_record = new HashMap<>();
        all_W = new HashMap<>();
        two_hops_neighbors = new HashSet[n];
        for (int i = 0; i < n; i++) {
            two_hops_neighbors[i] = new HashSet<>();
        }
    }

    /**
     * 计算所有顶点的two-hops邻居集合，用在后续的initialPhase计算顶点对之间的Saving
     */
    private void computeTwoHopsNeighbors(){

        for (int u = 0; u < n; u++) {
            Queue<Integer> queue = new LinkedList<>();
            int[] visited = new int[n];
            queue.offer(u);
            while (!queue.isEmpty()) {
                Integer v = queue.poll();
                if(visited[v] >= 3) break;
                two_hops_neighbors[u].add(v);
                int[] neighbors = Gr.successorArray(v);
                for (int i = 0; i < neighbors.length; i++) {
                    if (visited[neighbors[i]] == 0 && neighbors[i] != u) {
                        visited[neighbors[i]] = visited[v] + 1;
                        queue.offer(neighbors[i]);
                    }
                }
            }
        }
    }

    /**
     * 计算一个超点的two-hops超点邻居集合
     * @param super_node_id 超点的编号
     * @return
     */
    private Set<Integer> computeTwoHopsNeighbors(int super_node_id) {
        if(I[super_node_id] == -1) return new HashSet<>();
        int[] nodes = recoverSuperNode(super_node_id);
        Set<Integer> result = new HashSet<>();
        for (int u : nodes) {
            for (int v : two_hops_neighbors[u]) {
                if(result.contains(S[v])) continue;
                if(S[v] != -1) result.add(S[v]);
            }
        }
        return result;
    }

    /**
     * 某个顶点对进入优先队列，需要先判断是否已经在队列里面
     *
     * @param p
     */
    private void inQueueH(NodesPair p) {
        if (!H_record.containsKey(p.A))
            H_record.put(p.A, new HashSet<>());
        if (!H_record.containsKey(p.B))
            H_record.put(p.B, new HashSet<>());
        if (!H_record.get(p.A).contains(p.B)) {
            H_record.get(p.A).add(p.B);
            H_record.get(p.B).add(p.A);
        } else {
            H.remove(p);
        }
        H.offer(p);
    }

    /**
     * 某个顶点对出队，需要判断是否在队列里面
     *
     * @param p
     */
    private void outQueueH(NodesPair p) {
        if (H_record.get(p.A).contains(p.B)) {
            H.remove(p);
            H_record.get(p.A).remove(p.B);
            H_record.get(p.B).remove(p.A);
        }
    }


    private void processAffectedPairs(double threshold, NodesPair p) {
        HashSet<NodesPair> affected_pairs = new HashSet<>();

        // 处理优先队列H中包含超点p.A的所有顶点对
        if (H_record.containsKey(p.A)) {
            HashSet<Integer> contain_A = new HashSet<>(H_record.get(p.A));
            for (Integer C : contain_A) {
                NodesPair temp = new NodesPair(p.A, C);
                H.remove(temp);
                H_record.get(temp.A).remove(temp.B);
                H_record.get(temp.B).remove(temp.A);
                affected_pairs.add(temp);
            }
        }

        // 处理优先队列H中包含超点p.B的所有顶点对
        if (H_record.containsKey(p.B)) {
            HashSet<Integer> contain_B = new HashSet<>(H_record.get(p.B));
            for (Integer C : contain_B) {
                NodesPair temp_old = new NodesPair(p.B, C);
                H.remove(temp_old);
                H_record.get(temp_old.A).remove(temp_old.B);
                H_record.get(temp_old.B).remove(temp_old.A);
                NodesPair temp_new = new NodesPair(p.A, C);
                affected_pairs.add(temp_new);
            }
        }

        // 处理优先队列H中不包含p.A或p.B的顶点对，这些可能被影响到
        HashSet<Integer> affected_super_nodes = new HashSet<>();
        for (Integer u : all_W.get(p.A).keySet()) {
            affected_super_nodes.add(S[u]);
        }
        for (Integer W : affected_super_nodes) {
            for (Integer V : affected_super_nodes) {
                if (W == p.A || V == p.A || W == p.B || V == p.B) continue;
                if (V <= W) continue;
                NodesPair temp = new NodesPair(W, V);
                if (H_record.containsKey(temp.A) && H_record.get(temp.A).contains(temp.B)) {
                    H.remove(temp);
                    H_record.get(temp.A).remove(temp.B);
                    H_record.get(temp.B).remove(temp.A);
                }
                if(I[temp.A] == -1 || I[temp.B] == -1) continue;
                affected_pairs.add(temp);
            }
        }

        // 处理所有被影响到的顶点对
        for (NodesPair temp : affected_pairs) {
            temp.saving = computeSaving(all_W.get(temp.A), all_W.get(temp.B), temp.A, temp.B);
            if(temp.saving < threshold || temp.saving < 0.02 ) continue;
            inQueueH(temp);
        }
    }

    /**
     * 顶点初始化的阶段，Greedy算法需要进行重载
     */
    public double initialPhase(double threshold) {
        System.out.println("# Initial Phase");
        long startTime = System.currentTimeMillis();

        for (int A = 0; A < n; A++) {
            if(I[A] == -1) continue;
            Set<Integer> two_hops_supernode = computeTwoHopsNeighbors(A);
            for (Integer B : two_hops_supernode) {
                if(B <= A || I[B] == -1) continue;
                if((H_record.containsKey(A) && H_record.get(A).contains(B)) || (H_record.containsKey(B) && H_record.get(B).contains(A))) continue;
                NodesPair p = new NodesPair(A, B);
                if (!all_W.containsKey(p.A)) {
                    all_W.put(p.A, createW(p.A));
                }
                if (!all_W.containsKey(p.B)) {
                    all_W.put(p.B, createW(p.B));
                }
                if(all_W.get(p.A).size() == 0 || all_W.get(p.B).size() == 0) continue;
                p.saving = computeSaving(all_W.get(p.A), all_W.get(p.B), p.A, p.B);
                if(p.saving < threshold || p.saving < 0.02) continue;
                inQueueH(p);
//                if (threshold == 0.4) {
//                    System.out.println(p);
//                }
            }
        }

        return (System.currentTimeMillis() - startTime) / 1000.0;
    }

    /**
     * 顶点合并的阶段，SWeG和LDME算法需要进行重载
     *
     * @param threshold 合并阶段的阈值，低于阈值的顶点对不合并
     */
    public double mergePhase(double threshold) {
        System.out.println("# Merge Phase");
        System.out.println(String.format("Threshold=%5f", threshold));
        long startTime = System.currentTimeMillis();

        while (!H.isEmpty()) {
            // 队首出队，开始合并超点 p.A 和 p.B
            NodesPair p = H.poll();
            if(H_record.containsKey(p.A) && H_record.get(p.A).contains(p.B)) H_record.get(p.A).remove(p.B);
            if(H_record.containsKey(p.B) && H_record.get(p.B).contains(p.A)) H_record.get(p.B).remove(p.A);
            HashMap<Integer, Integer> w_update = updateW(all_W.get(p.A), all_W.get(p.B));
            all_W.replace(p.A, w_update);
            all_W.remove(p.B);
            updateSuperNode(p.A, p.B);
//            if (threshold == 0.0) {
//                System.out.println("Merge " + p.A + " and " + p.B + " with saving="+p.saving);
//            }
            // 处理被影响到的顶点对
            processAffectedPairs(threshold, p);
        }

        return (System.currentTimeMillis() - startTime) / 1000.0;
    }

    /**
     * @param iteration              迭代次数
     * @param print_iteration_offset 每执行多少次迭代就进行一次 encode 和 evaluate 进行结果输出
     */
    @Override
    public void run(int iteration, int print_iteration_offset){
        System.out.println("----------------------------------- Greedy ALGORITHM ----------------------------------------");
        long startTime = System.currentTimeMillis();
        computeTwoHopsNeighbors();
        System.out.println("Compute all two-hops neighbors takes " + ((System.currentTimeMillis()-startTime)/1000.0) + " seconds");
        for (int it = 1; it <= iteration; it++) {
            System.out.println("\n------------------------- ITERATION " + it);
//            double threshold = 1 / ((it + 1) * 1.0);
            double threshold = 0.5 - it * 0.05;
            System.out.println(String.format("@Time: %5f seconds", initialPhase(threshold)));
            System.out.println("After Initial phase, H.size():" + H.size());
            System.out.println(String.format("@Time: %5f seconds", mergePhase(threshold)));
            System.out.println("After Merge phase, H.size():" + H.size());
            if (it % print_iteration_offset == 0) {
//                System.out.println(String.format("@Time: %5f seconds", encodePhase()));
                System.out.println(String.format("@Time: %5f seconds", encodePhase_new()));
                evaluatePhase();
            }
        }
    }
}
