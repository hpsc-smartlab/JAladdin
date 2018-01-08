package it.uniparthenope;

import it.uniparthenope.Parser.MyJSONParser;

import java.io.Serializable;

public class DepartureParameters implements Serializable {
    private long year;
    private long month;
    private long day;
    private long hour;
    private long min;
    private long now_flag;

    public DepartureParameters(){
        //getting data from datetime_pars.json parsing
        MyJSONParser parser = new MyJSONParser("datetime_pars.json");
        this.year = parser.getValueAsLong("year");
        this.month = parser.getValueAsLong("month");
        this.day = parser.getValueAsLong("day");
        this.hour = parser.getValueAsLong("hour");
        this.min = parser.getValueAsLong("min");
        this.now_flag = 0;
    }

    public long getYear() {
        return year;
    }

    public long getMonth() {
        return month;
    }

    public long getDay() {
        return day;
    }

    public long getHour() {
        return hour;
    }

    public long getMin() {
        return min;
    }

    public long getNow_flag() {
        return now_flag;
    }

    public void setYear(long year) {
        this.year = year;
    }

    public void setMonth(long month) {
        this.month = month;
    }

    public void setDay(long day) {
        this.day = day;
    }

    public void setHour(long hour) {
        this.hour = hour;
    }

    public void setMin(long min) {
        this.min = min;
    }

    public void setNow_flag(long now_flag) {
        this.now_flag = now_flag;
    }
}
