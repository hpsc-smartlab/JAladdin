package it.uniparthenope.Boxing;

import java.util.ArrayList;

public class vessel_ResponseResults {
    private double[][] ship_v_LUT;
    private ArrayList<Double> H_array_m;
    private double polar_twa;
    private double polar_tws;

    public vessel_ResponseResults(double[][] ship_v_LUT, ArrayList<Double> H_array_m, double polar_twa, double polar_tws){
        this.ship_v_LUT = ship_v_LUT;
        this.H_array_m = H_array_m;
        this.polar_twa = polar_twa;
        this.polar_tws = polar_tws;
    }


    public double[][] getShip_v_LUT() {
        return ship_v_LUT;
    }



    public ArrayList<Double> getH_array_m() {
        return H_array_m;
    }

    public double getPolar_twa() {
        return polar_twa;
    }

    public double getPolar_tws() {
        return polar_tws;
    }
}
