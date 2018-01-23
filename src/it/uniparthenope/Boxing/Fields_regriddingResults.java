package it.uniparthenope.Boxing;


import it.uniparthenope.Parser.JSONManager;
import org.json.simple.parser.ParseException;

import java.io.IOException;

public class Fields_regriddingResults {
    private int[] time_steps;
    private double[][][] VTDH_Inset;
    private double[][][] VTPK_Inset;
    private double[][][] VDIR_Inset;
    private double[][][] windMAGN_Inset;
    private double[][][] windDIR_Inset;

    public Fields_regriddingResults(int[] time_steps, double[][][] VTDH_Inset, double[][][] VTPK_Inset, double[][][] VDIR_Inset,
                                    double[][][] windMAGN_Inset, double[][][] windDIR_Inset){
        this.time_steps = time_steps;
        this.VTDH_Inset = VTDH_Inset;
        this.VTPK_Inset = VTPK_Inset;
        this.VDIR_Inset = VDIR_Inset;
        this.windMAGN_Inset = windMAGN_Inset;
        this.windDIR_Inset = windDIR_Inset;
    }

    public Fields_regriddingResults(boolean flag) throws IOException, ParseException {
        if(flag==true){
            JSONManager reader = new JSONManager();
            reader.initReading("SerializedObjects/Fields_regridding.json");
            time_steps = reader.retrieveIntArray("time_steps");
            VTDH_Inset = reader.retrieveDouble3D("VTDH_Inset");
            VTPK_Inset = reader.retrieveDouble3D("VTPK_Inset");
            VDIR_Inset = reader.retrieveDouble3D("VDIR_Inset");
            windMAGN_Inset = reader.retrieveDouble3D("windMAGN_Inset");
            windDIR_Inset = reader.retrieveDouble3D("windDIR_Inset");
            reader.dispose();
        }
    }

    public void saveState() throws IOException {
        JSONManager writer = new JSONManager();
        writer.initWriting("SerializedObjects/Fields_regridding.json");
        writer.putIntArray("time_steps", time_steps);
        writer.putDouble3D("VTDH_Inset", VTDH_Inset);
        writer.putDouble3D("VTPK_Inset", VTPK_Inset);
        writer.putDouble3D("VDIR_Inset", VDIR_Inset);
        writer.putDouble3D("windMAGN_Inset", windMAGN_Inset);
        writer.putDouble3D("windDIR_Inset", windDIR_Inset);
        writer.dispose();
    }

    public int[] getTime_steps() {
        return time_steps;
    }

    public double[][][] getVTDH_Inset() {
        return VTDH_Inset;
    }

    public double[][][] getVTPK_Inset() {
        return VTPK_Inset;
    }

    public double[][][] getVDIR_Inset() {
        return VDIR_Inset;
    }

    public double[][][] getWindMAGN_Inset() {
        return windMAGN_Inset;
    }

    public double[][][] getWindDIR_Inset() {
        return windDIR_Inset;
    }
}
