package it.uniparthenope.Boxing;

import java.util.ArrayList;

public class nodes_free_form_barrierResults {
    private int Nfn;
    private ArrayList<Integer> free_nodes;
    private int red_fact;

    public nodes_free_form_barrierResults(int Nfn, ArrayList<Integer> free_nodes, int red_fact){
        this.Nfn = Nfn;
        this.free_nodes = free_nodes;
        this.red_fact = red_fact;
    }

    public int getNfn() {
        return Nfn;
    }

    public ArrayList<Integer> getFree_nodes() {
        return free_nodes;
    }

    public int getRed_fact() {
        return red_fact;
    }
}
