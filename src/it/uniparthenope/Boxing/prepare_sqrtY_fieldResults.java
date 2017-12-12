package it.uniparthenope.Boxing;

public class prepare_sqrtY_fieldResults {
    double[][][] VTDH_b;
    double[][][] VDIR_b;

    public prepare_sqrtY_fieldResults(double[][][] VTDH_b, double[][][] VDIR_b){
        this.VTDH_b = VTDH_b;
        this.VDIR_b = VDIR_b;
    }

    public double[][][] getVTDH_b() {
        return VTDH_b;
    }

    public double[][][] getVDIR_b() {
        return VDIR_b;
    }
}
