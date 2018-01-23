package it.uniparthenope;


import it.uniparthenope.Parser.JSONManager;
import org.json.simple.parser.ParseException;

import java.io.IOException;

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

    public EnvironmentalFields(boolean flag) throws IOException, ParseException {
        if(flag == true){
            JSONManager reader = new JSONManager();
            reader.initReading("SerializedObjects/EnvironmentalFields.json");
            wave = reader.retrieveLong("wave");
            current = reader.retrieveLong("current");
            stokes_drift = reader.retrieveLong("stokes_drift");
            wind = reader.retrieveLong("wind");
            deepWaterApprox = reader.retrieveLong("deepWaterApprox");
            analytic = reader.retrieveLong("analytic");
            reader.dispose();
        }
    }

    public void saveState() throws IOException {
        JSONManager writer = new JSONManager();
        writer.initWriting("SerializedObjects/EnvironmentalFields.json");
        writer.putLong("wave", wave);
        writer.putLong("current", current);
        writer.putLong("stokes_drift", stokes_drift);
        writer.putLong("wind", wind);
        writer.putLong("deepWaterApprox", deepWaterApprox);
        writer.putLong("analytic", analytic);
        writer.dispose();
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
