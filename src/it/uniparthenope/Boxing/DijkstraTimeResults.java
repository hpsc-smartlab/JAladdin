package it.uniparthenope.Boxing;

import java.util.LinkedList;

public class DijkstraTimeResults extends Dijkstra2DResults {
    double[] partial_times;

    public DijkstraTimeResults(double tdep_cost, LinkedList<Integer>  tdep_path, double[] tdep_partial_times){
        this.cost = tdep_cost;
        this.path = tdep_path;
        this.partial_times = tdep_partial_times;
    }

    public double[] getPartial_times() {
        return partial_times;
    }
}
