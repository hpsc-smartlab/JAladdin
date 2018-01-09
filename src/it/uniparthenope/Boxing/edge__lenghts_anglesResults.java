package it.uniparthenope.Boxing;

public class edge__lenghts_anglesResults {
    double[] theta_grid;
    double[] edge_lenght;
    double[] xm;
    double[] ym;

    public edge__lenghts_anglesResults(double[] theta_grid, double[] edge_lenght, double[] xm, double[] ym) {
        this.theta_grid = theta_grid;
        this.edge_lenght = edge_lenght;
        this.xm = xm;
        this.ym = ym;
    }

    public double[] getTheta_grid() {
        return theta_grid;
    }

    public double[] getEdge_lenght() {
        return edge_lenght;
    }

    public double[] getXm() {
        return xm;
    }

    public double[] getYm() {
        return ym;
    }
}
