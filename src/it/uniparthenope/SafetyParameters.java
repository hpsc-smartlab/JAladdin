package it.uniparthenope;

import it.uniparthenope.Parser.JSONManager;
import it.uniparthenope.Parser.MyJSONParser;
import org.json.simple.parser.ParseException;

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;

public class SafetyParameters implements Serializable {
    private long par_rolling;
    private long pureLossStab;
    private long surfRiding;
    private ArrayList<Long> criteria;

    public SafetyParameters(){
        //getting data from safety_pars.json parsing
        MyJSONParser parser = new MyJSONParser("safety_pars.json");
        this.par_rolling = parser.getValueAsLong("par_rolling");
        this.pureLossStab = parser.getValueAsLong("pureLossStab");
        this.surfRiding = parser.getValueAsLong("surfRiding");
        this.criteria = new ArrayList<Long>();
    }

    public SafetyParameters(boolean flag) throws IOException, ParseException{
        if(flag==true){
            JSONManager reader = new JSONManager();
            reader.initReading("SerializedObjects/SafetyParameters.json");
            par_rolling = reader.retrieveLong("par_rolling");
            pureLossStab = reader.retrieveLong("pureLossStab");
            surfRiding = reader.retrieveLong("surfRiding");
            criteria = reader.retrieveLongArrayList("criteria");
            reader.dispose();
        }
    }

    public void saveState() throws IOException {
        JSONManager writer = new JSONManager();
        writer.initWriting("SerializedObjects/SafetyParameters.json");
        writer.putLong("par_rolling", par_rolling);
        writer.putLong("pureLossStab", pureLossStab);
        writer.putLong("surfRiding", surfRiding);
        writer.putLongArrayList("criteria", criteria);
        writer.dispose();
    }

    public long getPar_rolling() {
        return par_rolling;
    }

    public long getPureLossStab() {
        return pureLossStab;
    }

    public long getSurfRiding() {
        return surfRiding;
    }

    public ArrayList<Long> getCriteria() {
        return criteria;
    }

    public void setCriteria(){
        this.criteria = new ArrayList<Long>();
        this.criteria.add(new Long(this.par_rolling));
        this.criteria.add(new Long(this.pureLossStab));
        this.criteria.add(new Long(this.surfRiding));
    }
}
