package it.uniparthenope;

import it.uniparthenope.Parser.JSONManager;
import org.json.simple.parser.ParseException;

import java.io.IOException;

public class Fstats {
    private double wheight_min;
    private double wheight_max;
    private double bathy_min;
    private double bathy_max;
    private double wperiod_min;
    private double wperiod_max;
    private double wheight_avg;
    private double wlength_min;
    private double wlenght_max;

    public Fstats(){}

    public Fstats(boolean flag) throws IOException, ParseException {
        if(flag==true){
            JSONManager reader = new JSONManager();
            reader.initReading("SerializedObjects/Fstats.json");
            wheight_min = reader.retrieveDouble("wheight_min");
            wheight_max = reader.retrieveDouble("wheight_max");
            bathy_min = reader.retrieveDouble("bathy_min");
            bathy_max = reader.retrieveDouble("bathy_max");
            wperiod_min = reader.retrieveDouble("wperiod_min");
            wperiod_max = reader.retrieveDouble("wperiod_max");
            wheight_avg = reader.retrieveDouble("wheight_avg");
            wlength_min = reader.retrieveDouble("wlength_min");
            wlenght_max = reader.retrieveDouble("wlenght_max");
            reader.dispose();
        }
    }

    public void saveState() throws IOException {
        JSONManager writer = new JSONManager();
        writer.initWriting("SerializedObjects/Fstats.json");
        writer.putDouble("wheight_min", wheight_min);
        writer.putDouble("wheight_max", wheight_max);
        writer.putDouble("bathy_min", bathy_min);
        writer.putDouble("wperiod_min", wperiod_min);
        writer.putDouble("wperiod_max", wperiod_max);
        writer.putDouble("wheight_avg", wheight_avg);
        writer.putDouble("wlength_min", wlength_min);
        writer.putDouble("wlenght_max", wlenght_max);
        writer.dispose();
    }

    public double getWheight_min() {
        return wheight_min;
    }

    public void setWheight_min(double wheight_min) {
        this.wheight_min = wheight_min;
    }

    public double getWheight_max() {
        return wheight_max;
    }

    public void setWheight_max(double wheight_max) {
        this.wheight_max = wheight_max;
    }

    public double getBathy_min() {
        return bathy_min;
    }

    public void setBathy_min(double bathy_min) {
        this.bathy_min = bathy_min;
    }

    public double getBathy_max() {
        return bathy_max;
    }

    public void setBathy_max(double bathy_max) {
        this.bathy_max = bathy_max;
    }

    public double getWheight_avg() {
        return wheight_avg;
    }

    public void setWheight_avg(double wheight_avg) {
        this.wheight_avg = wheight_avg;
    }

    public double getWperiod_min() {
        return wperiod_min;
    }

    public void setWperiod_min(double wperiod_min) {
        this.wperiod_min = wperiod_min;
    }

    public double getWperiod_max() {
        return wperiod_max;
    }

    public void setWperiod_max(double wperiod_max) {
        this.wperiod_max = wperiod_max;
    }

    public double getWlength_min() {
        return wlength_min;
    }

    public void setWlength_min(double wlength_min) {
        this.wlength_min = wlength_min;
    }

    public double getWlenght_max() {
        return wlenght_max;
    }

    public void setWlenght_max(double wlenght_max) {
        this.wlenght_max = wlenght_max;
    }
}
