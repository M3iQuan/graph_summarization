package graph_summarization;

/**
 * 用于记录组别的分布情况，即大小，数量和合并时间
 */
public class Record {
    Integer size;  // 组内顶点的数量
    Integer num;  // 组的数量
    Integer success;  // 成功合并的数量
    Double time;  // 组内合并时间

    public Record(Integer size, Integer num, Integer success, Double time){
        this.size = size;
        this.num = num;
        this.success = success;
        this.time = time;
    }

    public void add(Integer num, Integer success, Double time) {
        this.addNum(num);
        this.addSuccess(success);
        this.addTime(time);
    }

    public void addNum(Integer num) {
        this.num += num;
    }

    public void addSuccess(Integer success) {
        this.success += success;
    }

    public void addTime(Double time) {
        this.time += time;
    }

    public String toString() {
        return this.size + "," + this.num + "," + this.success + "," + String.format("%5f", this.time);
    }
}
