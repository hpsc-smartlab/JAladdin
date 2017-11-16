package it.uniparthenope;

public class DepartureParameters {
    private long year;
    private long month;
    private long day;
    private long hour;
    private long min;

    public DepartureParameters(){
        //getting data from datetime_pars.json parsing
        MyJSONParser parser = new MyJSONParser("datetime_pars.json");
        this.year = parser.getValueAsLong("year");
        this.month = parser.getValueAsLong("month");
        this.day = parser.getValueAsLong("day");
        this.hour = parser.getValueAsLong("hour");
        this.min = parser.getValueAsLong("min");
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
}
