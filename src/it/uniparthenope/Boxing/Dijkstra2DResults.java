package it.uniparthenope.Boxing;

import java.util.ArrayList;

public class Dijkstra2DResults {
    double cost;
    ArrayList<Integer> path;

    public Dijkstra2DResults(double cost, ArrayList<Integer> path){
        this.cost = cost;
        this.path = (ArrayList<Integer>) path.clone();
    }
}
