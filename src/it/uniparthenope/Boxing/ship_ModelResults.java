package it.uniparthenope.Boxing;

import java.util.ArrayList;

public class ship_ModelResults {
    private double[][] ship_v_LUT;
    private ArrayList<Double> H_array_m;

    public ship_ModelResults(double[][] ship_v_LUT, ArrayList<Double> H_array_m){
        this.ship_v_LUT = ship_v_LUT;
        this.H_array_m = H_array_m;
    }

    public double[][] getShip_v_LUT() {
        return ship_v_LUT;
    }

    public ArrayList<Double> getH_array_m() {
        return H_array_m;
    }
}
