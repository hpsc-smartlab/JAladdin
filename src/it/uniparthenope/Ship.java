package it.uniparthenope;

public class Ship {
    private int sailType;
    private double etaEngine;

    public Ship(){
        this.sailType = 4; //sailboats
        this.etaEngine = 0.7;

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
