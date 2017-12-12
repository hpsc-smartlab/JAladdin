package it.uniparthenope.Boxing;

public class fake_waveFieldsResults {
    int[] wave_origtimes;
    double[] lon_wave;
    double[] lat_wave;
    double[][][] VTDH;
    double[][][] VTPK;
    double[][][] VDIR;

    public fake_waveFieldsResults(int[] wave_origtimes, double[] lon_wave, double[] lat_wave, double[][][] VTDH, double[][][] VTPK, double[][][] VDIR){
        this.wave_origtimes = wave_origtimes;
        this.lon_wave = lon_wave;
        this.lat_wave = lat_wave;
        this.VTDH = VTDH;
        this.VTPK = VTPK;
        this.VDIR = VDIR;
    }

    public int[] getWave_origtimes() {
        return wave_origtimes;
    }

    public double[] getLon_wave() {
        return lon_wave;
    }

    public double[] getLat_wave() {
        return lat_wave;
    }

    public double[][][] getVTDH() {
        return VTDH;
    }

    public double[][][] getVTPK() {
        return VTPK;
    }

    public double[][][] getVDIR() {
        return VDIR;
    }
}
