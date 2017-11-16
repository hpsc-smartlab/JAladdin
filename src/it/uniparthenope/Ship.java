package it.uniparthenope;

public class Ship {
    /*defined in settings.m*/
    private int sailType;
    private double etaEngine;
    /***********************/
    private int vessType;
    private int sailClass;
    private int P_max_hp;
    private double maxv;
    private double length;
    private double beam;
    private double draught;
    private double obs_roll_period;

    public Ship(){
        this.sailType = 4; //sailboats
        this.etaEngine = 0.7;
    }

    public void LoadVesselParameters(){
        //getting data from ship_pars.json parsing
        MyJSONParser parser = new MyJSONParser("ship_pars.json");
        this.vessType = parser.getValueAsInt("vessType");
        this.sailClass = parser.getValueAsInt("sailClass");
        this.P_max_hp = parser.getValueAsInt("P_max_hp");
        this.maxv = parser.getValueAsDouble("maxv");
        this.length = parser.getValueAsDouble("length");
        this.beam = parser.getValueAsDouble("beam");
        this.draught = parser.getValueAsDouble("draught");
        this.obs_roll_period = parser.getValueAsDouble("obs_roll_period");
    }

    public int getVessType() {
        return vessType;
    }

    public int getSailClass() {
        return sailClass;
    }

    public int getP_max_hp() {
        return P_max_hp;
    }

    public double getMaxv() {
        return maxv;
    }

    public double getLength() {
        return length;
    }

    public double getBeam() {
        return beam;
    }

    public double getDraught() {
        return draught;
    }

    public double getObs_roll_period() {
        return obs_roll_period;
    }

    public double getEtaEngine() {
        return etaEngine;
    }

    public double getSailType(){
        return this.sailType;
    }

    public void setEtaEngine(double etaEngine) {
        this.etaEngine = etaEngine;
    }

    public void setSailType(int sailType){
        this.sailType = sailType;
    }
}
