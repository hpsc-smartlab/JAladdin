package it.uniparthenope.Boxing;

import java.util.ArrayList;

public class readout_bathyResults {
    ArrayList<Double> lat;
    ArrayList<Double> lon;
    Double[][] z;

    public readout_bathyResults(ArrayList<Double> lat, ArrayList<Double> lon, Double[][] z){
        this.lat = lat;
        this.lon = lon;
        this.z = z;
    }

    public ArrayList<Double> getLat() {
        return lat;
    }

    public ArrayList<Double> getLon() {
        return lon;
    }

    public Double[][] getZ() {
        return z;
    }
}
