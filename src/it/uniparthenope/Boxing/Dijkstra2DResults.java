package it.uniparthenope.Boxing;

import java.util.LinkedList;

public class Dijkstra2DResults {
    protected double cost;
    protected LinkedList<Integer> path;

    public Dijkstra2DResults(){};

    public Dijkstra2DResults(double cost, LinkedList<Integer>  path){
        this.cost = cost;
        this.path = path;
    }

    public double getCost() {
        return cost;
    }

    public LinkedList<Integer>  getPath() {
        return path;
    }
}
