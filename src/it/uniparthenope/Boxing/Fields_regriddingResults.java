package it.uniparthenope.Boxing;

import java.io.Serializable;

public class Fields_regriddingResults implements Serializable {
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
