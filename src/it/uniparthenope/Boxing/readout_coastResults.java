package it.uniparthenope.Boxing;

import java.util.ArrayList;

public class readout_coastResults {
    private ArrayList<Double> lon_ext;
    private ArrayList<Double> lat_ext;
    private ArrayList<Double> lon_int;
    private ArrayList<Double> lat_int;

    public readout_coastResults(ArrayList<Double> lon_ext, ArrayList<Double> lat_ext, ArrayList<Double> lon_int, ArrayList<Double> lat_int) {
        this.lon_ext = lon_ext;
        this.lat_ext = lat_ext;
        this.lon_int = lon_int;
        this.lat_int = lat_int;
    }

    public ArrayList<Double> getLon_ext() {
        return lon_ext;
    }

    public ArrayList<Double> getLat_ext() {
        return lat_ext;
    }

    public ArrayList<Double> getLon_int() {
        return lon_int;
    }

    public ArrayList<Double> getLat_int() {
        return lat_int;
    }
}
