package it.uniparthenope.Boxing;

import java.util.ArrayList;

public class grid_extreme_coordsResults {
    private ArrayList<Double> lat_red;
    private ArrayList<Double> lon_red;
    private Double[][] field_out;

    public grid_extreme_coordsResults(ArrayList<Double> lat_red, ArrayList<Double> lon_red, Double[][] field_out){
        this.lat_red = lat_red;
        this.lon_red = lon_red;
        this.field_out = field_out;
    }

    public ArrayList<Double> getLat_red() {
        return lat_red;
    }

    public ArrayList<Double> getLon_red() {
        return lon_red;
    }

    public Double[][] getField_out() {
        return field_out;
    }
}
