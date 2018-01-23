package it.uniparthenope.Boxing;

import it.uniparthenope.Parser.JSONManager;
import org.json.simple.parser.ParseException;

import java.io.IOException;
import java.util.ArrayList;

public class vessel_ResponseResults{
    private double[][] ship_v_LUT;
    private ArrayList<Double> H_array_m;
    private ArrayList<Double> polar_twa;
    private ArrayList<Double> polar_tws;

    public vessel_ResponseResults(boolean flag) throws IOException, ParseException{
        if(flag==true) {
            JSONManager reader = new JSONManager();
            reader.initReading("SerializedObjects/vesselResponse.json");
            this.ship_v_LUT = reader.retrieveDouble2D("ship_v_LUT");
            this.H_array_m = reader.retrieveDoubleArrayList("H_array_m");
            this.polar_twa = reader.retrieveDoubleArrayList("polar_twa");
            this.polar_tws = reader.retrieveDoubleArrayList("polar_tws");
            reader.dispose();
        }
    }

    public vessel_ResponseResults(double[][] ship_v_LUT, ArrayList<Double> H_array_m, ArrayList<Double> polar_twa, ArrayList<Double> polar_tws){
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

    public ArrayList<Double> getPolar_twa() {
        return polar_twa;
    }

    public ArrayList<Double> getPolar_tws() {
        return polar_tws;
    }

    public void saveState() throws IOException{
        JSONManager writer = new JSONManager();
        writer.initWriting("SerializedObjects/vesselResponse.json");
        writer.putDouble2D("ship_v_LUT", this.ship_v_LUT);
        writer.putDoubleArrayList("H_array_m", this.H_array_m);
        writer.putDoubleArrayList("polar_twa", this.polar_twa);
        writer.putDoubleArrayList("polar_tws", this.polar_tws);
        writer.dispose();
    }
}
