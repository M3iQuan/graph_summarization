package graph_summarization;

public class HashCodeAdjust {
    int codeSize;
    double[] code;

    public HashCodeAdjust(int codeSize) {
        this.codeSize = codeSize;
        code = new double[codeSize];
        for(int i=0; i < codeSize; i++) code[i] = 0.0;
    }

    public HashCodeAdjust(double[] code) {
        this.codeSize = code.length;
        this.code = code;
    }

    public boolean equals(HashCode otherCode) {
        for (int i = 0; i < codeSize; i++) {
            if (code[i] != otherCode.code[i]) {
                return false;
            }
        }
        return true;
    }

    public boolean unassigned(){
        for (int i = 0; i < codeSize; i++) {
            if(code[i] != 0){
                return false;
            }
        }
        return true;
    }

    public static int compare(HashCode a, HashCode b) {
        for (int i = 0; i < a.codeSize; i++) {
            if (a.code[i] != b.code[i]) {
                return a.code[i] - b.code[i];
            }
        }
        return 0;
    }

    public String toString(){
        StringBuilder stringBuilder = new StringBuilder();
        for (int i = 0; i < codeSize; i++) {
            stringBuilder.append(code[i]);
            stringBuilder.append(",");
        }
        return stringBuilder.toString();
    }
}
