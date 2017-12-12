package it.uniparthenope.Boxing;

public class fake_windFieldsResults {
    double[] wind_origTimes;
    double[] lat_wind;
    double[] lon_wind;
    double[][][] U10m;
    double[][][] V10m;

    public fake_windFieldsResults(double[] wind_origTimes, double[] lat_wind, double[] lon_wind, double[][][] U10m, double[][][] V10m){
        this.wind_origTimes = wind_origTimes;
        this.lat_wind = lat_wind;
        this.lon_wind = lon_wind;
        this.U10m = U10m;
        this.V10m = V10m;
    }

    public double[] getWind_origTimes() {
        return wind_origTimes;
    }

    public double[] getLat_wind() {
        return lat_wind;
    }

    public double[] getLon_wind() {
        return lon_wind;
    }

    public double[][][] getU10m() {
        return U10m;
    }

    public double[][][] getV10m() {
        return V10m;
    }
}
