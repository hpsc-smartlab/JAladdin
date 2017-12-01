package it.uniparthenope.Boxing;

public class mdata_gridResults {
    private double[][] xy;
    private double[][] xg;
    private double[][] yg;

    public mdata_gridResults(double[][] xy, double[][] xg, double[][] yg){
        this.xy = xy;
        this.xg = xg;
        this.yg = yg;
    }

    public double[][] getXy() {
        return xy;
    }

    public double[][] getXg() {
        return xg;
    }

    public double[][] getYg() {
        return yg;
    }
}
