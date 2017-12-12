package it.uniparthenope.Boxing;

public class readout_mWaveResults {
    private double[] lat;
    private double[] lon;
    private boolean wave_status;
    private String wave_filename;
    private double[][][] VTDH;
    private double[][][] VTPK;
    private double[][][] VDIR;

    public readout_mWaveResults(double[] lat, double[] lon, boolean wave_status, String wave_filename, double[][][] VTDH, double[][][] VTPK, double[][][] VDIR){
        this.lat = lat;
        this.lon = lon;
        this.wave_filename = wave_filename;
        this.wave_status = wave_status;
        this.VTDH = VTDH;
        this.VDIR = VDIR;
        this.VTPK = VTPK;
    }

    public double[] getLat() {
        return lat;
    }

    public double[] getLon() {
        return lon;
    }

    public boolean isWave_status() {
        return wave_status;
    }

    public String getWave_filename() {
        return wave_filename;
    }

    public double[][][] getVTHD() {
        return VTDH;
    }

    public double[][][] getVTPK() {
        return VTPK;
    }

    public double[][][] getVDIR() {
        return VDIR;
    }
}
