package it.uniparthenope.Boxing;

public class readout_envFieldsResults {
    double[] lat_wave;
    double[] lon_wave;
    double[] ecmwf_lat_wind;
    double[] ecmwf_lon_wind;
    double[] cosmo_lat_wind;
    double[] cosmo_lon_wind;
    double[][][] VTDH;
    double[][][] VTPK;
    double[][][] VDIR;
    double[][][] ecmwf_U10m;
    double[][][] ecmwf_V10m;
    double[][][] cosmo_U10m;
    double[][][] cosmo_V10m;
    int[] wave_origtimes;
    double[] ecmwf_wind_origTimes;
    double[] cosmo_wind_origTimes;

    public readout_envFieldsResults(double[] lat_wave, double[] lon_wave, double[] ecmwf_lat_wind, double[] ecmwf_lon_wind,
                                    double[] cosmo_lat_wind, double[] cosmo_lon_wind, double[][][] VTDH, double[][][] VTPK,
                                    double[][][] VDIR, double[][][] ecmwf_U10m, double[][][] ecmwf_V10m, double[][][] cosmo_U10m,
                                    double[][][] cosmo_V10m, int[] wave_origtimes, double[] ecmwf_wind_origTimes, double[] cosmo_wind_origTimes){

        this.lat_wave = lat_wave;
        this.lon_wave = lon_wave;
        this.ecmwf_lat_wind = ecmwf_lat_wind;
        this.ecmwf_lon_wind = ecmwf_lon_wind;
        this.cosmo_lat_wind = cosmo_lat_wind;
        this.cosmo_lon_wind = cosmo_lon_wind;
        this.VTDH = VTDH;
        this.VTPK = VTPK;
        this.VDIR = VDIR;
        this.ecmwf_U10m = ecmwf_U10m;
        this.ecmwf_V10m = ecmwf_V10m;
        this.cosmo_U10m = cosmo_U10m;
        this.cosmo_V10m = cosmo_V10m;
        this.wave_origtimes = wave_origtimes;
        this.ecmwf_wind_origTimes = ecmwf_wind_origTimes;
        this.cosmo_wind_origTimes = cosmo_wind_origTimes;
    }

    public double[] getLat_wave() {
        return lat_wave;
    }

    public double[] getLon_wave() {
        return lon_wave;
    }

    public double[] getEcmwf_lat_wind() {
        return ecmwf_lat_wind;
    }

    public double[] getEcmwf_lon_wind() {
        return ecmwf_lon_wind;
    }

    public double[] getCosmo_lat_wind() {
        return cosmo_lat_wind;
    }

    public double[] getCosmo_lon_wind() {
        return cosmo_lon_wind;
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

    public double[][][] getEcmwf_U10m() {
        return ecmwf_U10m;
    }

    public double[][][] getEcmwf_V10m() {
        return ecmwf_V10m;
    }

    public double[][][] getCosmo_U10m() {
        return cosmo_U10m;
    }

    public double[][][] getCosmo_V10m() {
        return cosmo_V10m;
    }

    public int[] getWave_origtimes() {
        return wave_origtimes;
    }

    public double[] getEcmwf_wind_origTimes() {
        return ecmwf_wind_origTimes;
    }

    public double[] getCosmo_wind_origTimes() {
        return cosmo_wind_origTimes;
    }
}
