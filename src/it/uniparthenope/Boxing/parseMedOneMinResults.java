package it.uniparthenope.Boxing;

import java.util.ArrayList;

public class parseMedOneMinResults {
    ArrayList<Double> lat;
    ArrayList<Double> lon;
    double[][] depth;

    public parseMedOneMinResults(ArrayList<Double> lat, ArrayList<Double> lon, double[][] depth){
        this.lat = lat;
        this.lon = lon;
        this.depth = depth;
    }

    public ArrayList<Double> getLat() {
        return lat;
    }

    public ArrayList<Double> getLon() {
        return lon;
    }

    public double[][] getDepth() {
        return depth;
    }
}
