package it.uniparthenope;

public class EnvironmentalFields {
    private long wave;
    private long current;
    private long stokes_drift;
    private long wind;
    private long deepWaterApprox;
    private long analytic;

    public EnvironmentalFields(long shipVessType, long shipSailType, long windModel){
        if(shipVessType != shipSailType){ //motorboat
            this.wave = 1;
            this.current = 0;
            this.stokes_drift = 0;
            this.wind = 0;
        } else { //sailboat
            this.wave = 0;
            this.current = 0;
            this.stokes_drift = 0;
            this.wind = windModel;
        }
        this.deepWaterApprox = 1; //1: deep water approx;  0:Fenton-McKee wave dispersion in T_E and lambda
        //Just for GMD paper s route #2 (1589_c) uncomment following line:
        //this.deepWaterApprox = 0;
        this.analytic = 0; //if analytic = 1, it is mandatory to set timedep_flag=0  and vessType=1,2,3 (motorboat) in ship_pars.json
    }

    public long getWave() {
        return wave;
    }

    public long getCurrent() {
        return current;
    }

    public long getStokes_drift() {
        return stokes_drift;
    }

    public long getWind() {
        return wind;
    }

    public long getDeepWaterApprox() {
        return deepWaterApprox;
    }

    public long getAnalytic() {
        return analytic;
    }
}
