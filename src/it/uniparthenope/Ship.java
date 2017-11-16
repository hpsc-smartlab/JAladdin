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
    private double roll_period;
    private double bottomDraughtfactor;
    private double TEsaturFactor;
    private double Rspectrum;
    private double Nvel;

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
        if( this.vessType != this.sailType){
            this.roll_period = this.obs_roll_period;
        } else {
            this.roll_period = 1;
        }
        this.bottomDraughtfactor = 1; //bottomDraughtfactor*draught is threshold for effect of sea bottom on enhanced frictional
        this.TEsaturFactor = 10; //max allowed TE in plots is ship.roll_period*ship.TEsaturFactor
        this.Rspectrum = 0;
        this.Nvel = 1; //Default value for steps in power reduction.
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

    public void setStepsInPowerReduction(long flag){
        if(flag==1){
            this.Nvel = 7;
        } else { this.Nvel = 1; }
        if (this.vessType == this.sailType){ //sailboat
            this.Nvel = 1;
        }
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

    public double getRoll_period() {
        return roll_period;
    }

    public double getBottomDraughtfactor() {
        return bottomDraughtfactor;
    }

    public double getTEsaturFactor() {
        return TEsaturFactor;
    }

    public double getRspectrum() {
        return Rspectrum;
    }

    public double getNvel() {
        return Nvel;
    }

    public void setVessType(long vessType) {
        this.vessType = vessType;
    }
}
