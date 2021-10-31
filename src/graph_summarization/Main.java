package graph_summarization;

import java.util.HashMap;

public class Main {

    public static void testSWeG(String basename, int iteration, int print_iteration_offset) throws Exception{
        // 使用范型的方式，声明一个父类Summary,指向一个SWeG算法对象
        Summary S = new SWeG(basename);
        // 调用run方法运行整个压缩算法
        S.run(iteration, print_iteration_offset);
    }

    public static void testLDME(String basename, int iteration, int print_iteration_offset, int signatureLength) throws Exception{
        // 使用范型的方式，声明一个父类Summary,指向一个SWeG算法对象
        Summary S = new LDME(basename, signatureLength);
        // 调用run方法运行整个压缩算法
        S.run(iteration, print_iteration_offset);
    }

    public static void testGreedy(String basename, int iteration, int print_iteration_offset) throws Exception{
        // 使用范型的方式，声明一个父类Summary,指向一个SWeG算法对象
        Summary S = new Greedy(basename);
        // 调用run方法运行整个压缩算法
        S.run(iteration, print_iteration_offset);
    }

    public static void main(String[] args) throws Exception{
        // 参数读取,一共有四个 basename iteration print_iteration_offset k(只有LDME算法有)
        String basename = args[0];
        int iteration = Integer.parseInt(args[1]);
        int print_iteration_offset = Integer.parseInt(args[2]);
        int signatureLength = Integer.parseInt(args[3]);

//        testGreedy(basename, iteration, print_iteration_offset);
        testLDME(basename, iteration, print_iteration_offset, signatureLength);
//        testSWeG(basename, iteration, print_iteration_offset);
    }
}
