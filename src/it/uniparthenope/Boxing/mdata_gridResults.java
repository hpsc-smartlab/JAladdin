package it.uniparthenope.Boxing;

public class mdata_gridResults {
    private Double[][] xy;
    private Double[][] xg;
    private Double[][] yg;

    public mdata_gridResults(Double[][] xy, Double[][] xg, Double[][] yg){
        this.xy = xy;
        this.xg = xg;
        this.yg = yg;
    }

    public Double[][] getXy() {
        return xy;
    }

    public Double[][] getXg() {
        return xg;
    }

    public Double[][] getYg() {
        return yg;
    }
}
