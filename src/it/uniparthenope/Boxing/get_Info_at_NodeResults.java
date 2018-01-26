package it.uniparthenope.Boxing;

public class get_Info_at_NodeResults {
    private double[] Dr_cum;
    private double[] v_opt;
    private double[] edge_delay;
    private double[] theta_opt;
    private double[] theta_VMC;

    public get_Info_at_NodeResults(double[] dr_cum, double[] v_opt, double[] edge_delay, double[] theta_opt, double[] theta_VMC) {
        Dr_cum = dr_cum;
        this.v_opt = v_opt;
        this.edge_delay = edge_delay;
        this.theta_opt = theta_opt;
        this.theta_VMC = theta_VMC;
    }

    public double[] getDr_cum() {
        return Dr_cum;
    }

    public double[] getV_opt() {
        return v_opt;
    }

    public double[] getEdge_delay() {
        return edge_delay;
    }

    public double[] getTheta_opt() {
        return theta_opt;
    }

    public double[] getTheta_VMC() {
        return theta_VMC;
    }
}
