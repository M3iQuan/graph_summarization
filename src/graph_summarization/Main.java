package graph_summarization;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.HashMap;
import java.util.Map;

public class Main {

    public static void testSWeG(String basename, int iteration, int print_iteration_offset, int signatureLength, int max_group_size, int hierarchical_k) throws Exception{
        // 使用范型的方式，声明一个父类Summary,指向一个SWeG算法对象
        Summary S = new SWeG(basename);
//        Summary S = new SWeG(basename, signatureLength);
        // 调用run方法运行整个压缩算法
        S.run(iteration, print_iteration_offset, max_group_size, hierarchical_k);
//        S.run(iteration, print_iteration_offset);
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

    public static Map<String, Object> processArgs(String[] args) {
        Map<String, Object> result = new HashMap<>();
        String method="LDME", basename="./data/cnr-2000/cnr-2000-sym";
        int iteration=60, print_iteration_offset=10, signature_length=0, max_group_size=0, hierarchical_k=0;
        for (String arg : args) {
            String[] arr = arg.split("=");
            String key = arr[0];
            switch (key) {
                case "dataset":
                    basename = arr[1];
                    break;
                case "method":
                    method = arr[1];
                    break;
                case "iteration":
                    iteration = Integer.parseInt(arr[1]);
                    break;
                case "print_iteration_offset":
                    print_iteration_offset = Integer.parseInt(arr[1]);
                    break;
                case "signature":
                    signature_length = Integer.parseInt(arr[1]);
                    break;
                case "max_group_size":
                    max_group_size = Integer.parseInt(arr[1]);
                    break;
                case "hierarchical_k":
                    hierarchical_k = Integer.parseInt(arr[1]);
                    break;
            }
        }
        result.put("method", method);
        result.put("dataset", basename);
        result.put("iteration", iteration);
        result.put("print_iteration_offset", print_iteration_offset);
        if (signature_length != 0) {
            result.put("signature_length", signature_length);
        }
        if (max_group_size != 0) {
            result.put("max_group_size", max_group_size);
        }
        if (hierarchical_k != 0) {
            result.put("hierarchical_k", hierarchical_k);
        }
        return result;
    }

    public static void main(String[] args) throws Exception{
        Summary S;

        // 日志对象，用于记录参数的读取和调用的方法
        Logger logger = LoggerFactory.getLogger(Main.class);

        // 参数读取, 参数默认设置
        Map<String, Object> arguments = processArgs(args);
        StringBuilder message = new StringBuilder("调用%s算法");
        if(!arguments.containsKey("max_group_size")){
            message.append(", 不限制小组顶点数量");
        } else{
            message.append(", 限制小组顶点数量");
        }
        if(!arguments.containsKey("hierarchical_k")){
            message.append(", 不进行层次分组");
        } else{
            message.append(", 进行层次分组");
        }
        String method = (String)arguments.get("method");
        String basename = (String) arguments.get("dataset");
        logger.info(String.format(message.toString(), method));
        message = new StringBuilder("程序运行参数: ");
        switch (method) {
            case "LDME":
                S = new LDME(basename, (Integer)arguments.get("signature_length"));
                break;
            case "SWeG":
                S = new SWeG(basename);
//                S = new SWeG(basename, (Integer)arguments.get("signature_length"));
                break;
            case "Greedy":
//                S = new Greedy(basename);
                break;
            default:
                throw new IllegalStateException("不存在这个方法: " + method);
        }
        for (String key : arguments.keySet()) {
           message.append(key).append("=").append(arguments.get(key)).append(" ");
        }
        logger.info(message.toString());
//        S.run(iteration, print_iteration_offset, max_group_size, hierarchical_k);
    }
}
