package it.uniparthenope;

import java.io.Serializable;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.TimeZone;

public class TemporalGrid implements Serializable {
    private double Dt;
    private long minNt;
    private double Nt;
    private String latest_date;
    private long depDateTime;
    private long wave_dep_TS;
    private long wind_dep_TS;
    private double maxNt;



    public TemporalGrid(){
        this.Dt = 1.00;//h
        this.minNt = 24;
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

    public void setMaxNt(double maxNt) {
        this.maxNt = maxNt;
    }
}
