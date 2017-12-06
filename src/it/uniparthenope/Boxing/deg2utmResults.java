package it.uniparthenope.Boxing;

public class deg2utmResults {
    private double[] x;
    private double[] y;
    private String[][] utmzone;

    public deg2utmResults(double[] x, double[] y, String[][] utmzone){
        this.x = x;
        this.y = y;
        this.utmzone = utmzone;
    }

    public double[] getX() {
        return x;
    }

    public double[] getY() {
        return y;
    }

    public String[][] getUtmzone() {
        return utmzone;
    }
}
