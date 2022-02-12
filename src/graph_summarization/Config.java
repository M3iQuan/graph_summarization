package graph_summarization;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * 程序运行的配置类
 */
public class Config {
    // 运行的数据集
    String dataset;
    // 运行的方法, 目前包含 SWeG 和 LDME
    String method;
    // 迭代次数
    int iteration;
    // 多少轮迭代后输出一次性能
    int print_iteration_offset;
    // 签名长度, 用于Local Sensitive Hash函数
    int signature_length;
    // 小组的最大顶点数量
    int max_group_size;
    // 是否采用层次划分, 如果为false采用顺序切割, 为true则采用层次切割
    boolean hierarchical;
    // 如果不为0, 则采用的是n*k
    int hierarchical_k;

    // 日志文件
    Logger logger_ = LoggerFactory.getLogger(getClass());

    /**
     * 默认构造函数
     */
    public Config(){

    }

    /**
     * 初始化参数
     */
    public void init(){
        logger_.info("没有输入参数, 程序采用默认的参数执行");
        method = "SWeG";
        dataset = "./data/cnr-2000/cnr-2000-sym";
        iteration = 60;
        print_iteration_offset = 10;
        signature_length = 0;
        max_group_size = 0;
        hierarchical = false;
        hierarchical_k = 0;
    }

    public void outputMessage(){
        String message = "程序运行参数: method=%s dataset=%s iteration=%d print_iteration_offset=%d ";
        logger_.info("程序运行参数");
    }

    /**
     * 处理程序输入的参数
     * @param args
     */
    public void parseArguments(String[] args){
        int sum = 0;
        for(String arg : args){
            String key = arg.split("=")[0];
            String value = arg.split("=")[1];
            switch (key){
                case "method":
                    method = value;
                    sum++;
                    break;
                case "dataset":
                    dataset = value;
                    sum++;
                    break;
                case  "iteration":
                    iteration = Integer.parseInt(value);
                    sum++;
                    break;
                case "print_iteration_offset":
                    print_iteration_offset = Integer.parseInt(value);
                    sum++;
                    break;
            }
        }


        if(args.length == 0){
            init();
        }
    }
}
