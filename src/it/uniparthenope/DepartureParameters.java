package it.uniparthenope;

import it.uniparthenope.Parser.JSONManager;
import it.uniparthenope.Parser.MyJSONParser;
import org.json.simple.parser.ParseException;

import java.io.IOException;


public class DepartureParameters {
    private long year;
    private long month;
    private long day;
    private long hour;
    private long min;
    private long now_flag;
    private String depDateTime;

    public DepartureParameters(boolean flag) throws IOException, ParseException{
        if(flag == true){
            JSONManager reader = new JSONManager();
            reader.initReading("SerializedObjects/DepartureParameters.json");
            year = reader.retrieveLong("year");
            month = reader.retrieveLong("month");
            day = reader.retrieveLong("day");
            hour = reader.retrieveLong("hour");
            min = reader.retrieveLong("min");
            now_flag = reader.retrieveLong("now_flag");
            depDateTime = reader.retrieveString("depDateTime");
            reader.dispose();
        }
    }

    public DepartureParameters(){
        //getting data from datetime_pars.json parsing
        MyJSONParser parser = new MyJSONParser("datetime_pars.json");
        this.year = parser.getValueAsLong("year");
        this.month = parser.getValueAsLong("month");
        this.day = parser.getValueAsLong("day");
        this.hour = parser.getValueAsLong("hour");
        this.min = parser.getValueAsLong("min");
        this.now_flag = 0;
        setDepDateTime();
    }

    public void saveState() throws IOException {
        JSONManager writer = new JSONManager();
        writer.initWriting("SerializedObjects/DepartureParameters.json");
        writer.putLong("year", year);
        writer.putLong("month", month);
        writer.putLong("day", day);
        writer.putLong("hour", hour);
        writer.putLong("min", min);
        writer.putLong("now_flag", now_flag);
        writer.putString("depDateTime", depDateTime);
        writer.dispose();
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

    public String getDepDateTime() {
        return depDateTime;
    }

    private void setDepDateTime(){
        String stringYear = ""+year;
        if((year-100)<100)//check year format (2015 or 15)
            stringYear="20"+year;
        String stringMonth = ""+month;
        if(month<10)
            stringMonth = "0"+month;
        String stringday = ""+day;
        if(day<10)
            stringday = "0"+day;
        String stringhour=""+hour;
        if(hour<10)
            stringhour = "0"+hour;
        String stringMin=""+min;
        if(min<10)
            stringMin="0"+min;

        this.depDateTime = stringYear+stringMonth+stringday+stringhour+stringMin;

    }

    public void setNow_flag(long now_flag) {
        this.now_flag = now_flag;
    }
}
