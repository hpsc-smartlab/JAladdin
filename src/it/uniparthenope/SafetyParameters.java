package it.uniparthenope;

import it.uniparthenope.Parser.MyJSONParser;

import java.util.ArrayList;

public class SafetyParameters {
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
