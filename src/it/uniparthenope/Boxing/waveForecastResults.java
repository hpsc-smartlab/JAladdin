package it.uniparthenope.Boxing;

public class waveForecastResults {
    int[] time;
    double[] latitude;
    double[] longitude;
    double[][][] VDIR;
    double[][][] VTDH;
    double[][][] VTPK;

    public waveForecastResults(double[][][] VDIR, double[][][] VTDH, double[][][] VTPK, int[] time, double[] latitude, double[] longitude){
        this.time = time;
        this.latitude = latitude;
        this.longitude = longitude;
        this.VDIR = VDIR;
        this.VTDH = VTDH;
        this.VTPK = VTPK;
    }

    public int[] getTime() {
        return time;
    }

    public double[] getLatitude() {
        return latitude;
    }

    public double[] getLongitude() {
        return longitude;
    }

    public double[][][] getVDIR() {
        return VDIR;
    }

    public double[][][] getVTDH() {
        return VTDH;
    }

    public double[][][] getVTPK() {
        return VTPK;
    }
}
