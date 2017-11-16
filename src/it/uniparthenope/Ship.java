package it.uniparthenope;

public class Ship {
    /*defined in settings.m*/
    private long sailType;
    private double etaEngine;
    /***********************/
    private long vessType;
    private long sailClass;
    private long P_max_hp;
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
        this.vessType = parser.getValueAsLong("vessType");
        this.sailClass = parser.getValueAsLong("sailClass");
        this.P_max_hp = parser.getValueAsLong("P_max_hp");
        this.maxv = parser.getValueAsDouble("maxv");
        this.length = parser.getValueAsDouble("length(m)");
        this.beam = parser.getValueAsDouble("beam(m)");
        this.draught = parser.getValueAsDouble("draught(m)");
        this.obs_roll_period = parser.getValueAsDouble("obs_roll_period(s)");
        if(!IntegrityCheck()){
            System.out.println("Ship integrity check violated.");
            System.exit(-1);
        }
    }

    private boolean IntegrityCheck(){
        boolean check = true;
        if(this.vessType != this.sailType){
            if( this.P_max_hp <= 0 ||
                    this.maxv <= 0 ||
                    this.length <= 0 ||
                    this.beam <=0 ||
                    this.draught <= 0 ||
                    this.obs_roll_period <= 0){
                check = false;
            }
        }
        return check;
    }

    public long getVessType() {
        return vessType;
    }

    public long getSailClass() {
        return sailClass;
    }

    public long getP_max_hp() {
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

    public long getSailType(){
        return this.sailType;
    }

}
