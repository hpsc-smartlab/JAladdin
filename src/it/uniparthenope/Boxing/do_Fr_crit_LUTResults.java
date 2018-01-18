package it.uniparthenope.Boxing;

public class do_Fr_crit_LUTResults {
    double[] steep_var;
    double[] Fr_LUT;

    public do_Fr_crit_LUTResults(double[] xvar, double[] fr_LUT) {
        this.steep_var = xvar;
        Fr_LUT = fr_LUT;
    }

    public double[] getSteep_var() {
        return steep_var;
    }

    public double[] getFr_LUT() {
        return Fr_LUT;
    }
}
