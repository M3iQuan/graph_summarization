package graph_summarization;

public class NodesPair implements Comparable<NodesPair> {
    int A;
    int B;
    double saving = 0;
    double jac_sim = 0;

    public NodesPair(int super_node_a, int super_node_b){
        if (super_node_a < super_node_b) {
            A = super_node_a;
            B = super_node_b;
        } else {
            A = super_node_b;
            B = super_node_a;
        }
    }

    public NodesPair(int super_node_a, int super_node_b, double s) {
        if (super_node_a < super_node_b) {
            A = super_node_a;
            B = super_node_b;
        } else {
            A = super_node_b;
            B = super_node_a;
        }
        saving = s;
    }

    @Override
    public int compareTo(NodesPair o) {
        if (o.saving == this.saving) {
            if (o.A == this.A) {
                return this.B - o.B;
            } else {
                return this.A - o.A;
            }
        } else if (o.saving > this.saving){
            return 1;
        }
        return -1;
    }

    @Override
    public boolean equals(Object obj) {
        NodesPair o = (NodesPair) obj;
        return ((this.A == o.A) && (this.B == o.B)) || ((this.A == o.B) && (this.B == o.A));
    }

    @Override
    public String toString(){
        return "A:" + this.A + ", B:" + this.B + ", saving:" + this.saving;
    }

}
