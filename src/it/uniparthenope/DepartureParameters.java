package it.uniparthenope;

public class DepartureParameters {
    private int year;
    private int month;
    private int day;
    private int hour;
    private int min;

    public DepartureParameters(){
        //getting data from datetime_pars.json parsing
        MyJSONParser parser = new MyJSONParser("datetime_pars.json");
        this.year = parser.getValueAsInt("year");
        this.month = parser.getValueAsInt("month");
        this.day = parser.getValueAsInt("day");
        this.hour = parser.getValueAsInt("hour");
        this.min = parser.getValueAsInt("min");
    }

    public int getYear() {
        return year;
    }

    public int getMonth() {
        return month;
    }

    public int getDay() {
        return day;
    }

    public int getHour() {
        return hour;
    }

    public int getMin() {
        return min;
    }
}
