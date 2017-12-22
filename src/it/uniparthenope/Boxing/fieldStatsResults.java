package it.uniparthenope.Boxing;

public class fieldStatsResults {
    double ecmwf_dir_avg;
    double ecmwf_dir_std;
    double cosmo_dir_avg;
    double cosmo_dir_std;

    public fieldStatsResults(double ecmwf_dir_avg, double ecmwf_dir_std, double cosmo_dir_avg, double cosmo_dir_std){
        this.cosmo_dir_avg = cosmo_dir_avg;
        this.cosmo_dir_std = cosmo_dir_std;
        this.ecmwf_dir_avg = ecmwf_dir_avg;
        this.ecmwf_dir_std = ecmwf_dir_std;
    }

    public double getEcmwf_dir_avg() {
        return ecmwf_dir_avg;
    }

    public double getEcmwf_dir_std() {
        return ecmwf_dir_std;
    }

    public double getCosmo_dir_avg() {
        return cosmo_dir_avg;
    }

    public double getCosmo_dir_std() {
        return cosmo_dir_std;
    }
}
