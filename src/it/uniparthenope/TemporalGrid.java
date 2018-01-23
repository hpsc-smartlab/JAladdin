package it.uniparthenope;


import it.uniparthenope.Parser.JSONManager;
import org.json.simple.parser.ParseException;

import java.io.IOException;

public class TemporalGrid{
    private double Dt;
    private long minNt;
    private double Nt;
    private String latest_date;
    private long depDateTime;
    private long wave_dep_TS;
    private long wind_dep_TS;
    private double maxNt;
    private double[] dayHrs;



    public TemporalGrid(){
        this.Dt = 1.00;//h
        this.minNt = 24;
    }

    public TemporalGrid(boolean flag) throws IOException, ParseException{
        if(flag==true){
            JSONManager reader = new JSONManager();
            reader.initReading("SerializedObjects/TemporalGrid.json");
            Dt = reader.retrieveDouble("Dt");
            minNt = reader.retrieveLong("minNt");
            Nt = reader.retrieveDouble("Nt");
            latest_date = reader.retrieveString("latest_date");
            depDateTime = reader.retrieveLong("depDateTime");
            wave_dep_TS = reader.retrieveLong("wave_dep_TS");
            wind_dep_TS = reader.retrieveLong("wind_dep_TS");
            maxNt = reader.retrieveDouble("maxNt");
            dayHrs = reader.retrieveDoubleArray("dayHrs");
            reader.dispose();
        }
    }

    public void saveState() throws IOException {
        JSONManager writer = new JSONManager();
        writer.initWriting("SerializedObjects/TemporalGrid.json");
        writer.putDouble("Dt", Dt);
        writer.putLong("minNt", minNt);
        writer.putDouble("Nt", Nt);
        writer.putString("latest_date", latest_date);
        writer.putLong("depDateTime", depDateTime);
        writer.putLong("wave_dep_TS", wave_dep_TS);
        writer.putLong("wind_dep_TS", wind_dep_TS);
        writer.putDouble("maxNt", maxNt);
        writer.putDoubleArray("dayHrs", dayHrs);
        writer.dispose();
    }

    public double getDt() {
        return Dt;
    }

    public long getMinNt() {
        return minNt;
    }

    public String getLatest_date() {
        return latest_date;
    }

    public void setLatest_date(String latest_date) {
        this.latest_date = latest_date;
    }

    public void setDepDateTime(long depDateTime){
        this.depDateTime = depDateTime;
    }

    public long getDepDateTime() {
        return depDateTime;
    }

    public long getWave_dep_TS() {
        return wave_dep_TS;
    }

    public void setWave_dep_TS(long wave_dep_TS) {
        this.wave_dep_TS = wave_dep_TS;
    }

    public long getWind_dep_TS() {
        return wind_dep_TS;
    }

    public void setWind_dep_TS(long wind_dep_TS) {
        this.wind_dep_TS = wind_dep_TS;
    }

    public void setNt(double nt) {
        Nt = nt;
    }

    public double getNt() {
        return Nt;
    }

    public double getMaxNt() {
        return maxNt;
    }

    public double[] getDayHrs() {
        return dayHrs;
    }

    public void setDayHrs(double[] dayHrs) {
        this.dayHrs = dayHrs;
    }

    public void setMaxNt(double maxNt) {
        this.maxNt = maxNt;
    }
}
