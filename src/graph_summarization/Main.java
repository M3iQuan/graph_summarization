package graph_summarization;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class Main {

    public static void testSWeG(String basename, int iteration, int print_iteration_offset, int signatureLength) throws Exception{
        // 使用范型的方式，声明一个父类Summary,指向一个SWeG算法对象
//        Summary S = new SWeG(basename);
        Summary S = new SWeG(basename, signatureLength);
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
        // 日志对象，用于记录参数的读取和调用的方法
        Logger logger = LoggerFactory.getLogger(Main.class);
        // 参数读取,一共有五个 method, basename iteration print_iteration_offset k(只有LDME算法有)
        String method = args[0];
        // basename是数据集的名字
        String basename = args[1];
        // iteration是迭代次数
        int iteration = Integer.parseInt(args[2]);
        // print_iteration_offset是每迭代多少次就进行一次EdgeEncode，打印一次压缩率
        int print_iteration_offset = Integer.parseInt(args[3]);
        // signatureLength是LSH的签名长度，即论文里的k
        int signatureLength = Integer.parseInt(args[4]);

        logger.info(String.format("调用%s算法, basename=%s, iteration=%d, print_iteration_offset=%d, signatureLength=%d", method, basename, iteration, print_iteration_offset, signatureLength));
        switch (method) {
            case "LDME":
                testLDME(basename, iteration, print_iteration_offset, signatureLength);
                break;
            case "SWeG":
                testSWeG(basename, iteration, print_iteration_offset, signatureLength);
                break;
            case "Greedy":
                testGreedy(basename, iteration, print_iteration_offset);
                break;
        }
    }
}
