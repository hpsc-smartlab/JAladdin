package it.uniparthenope;

public class Const {
    private double deg2rad;
    private double m2ft;
    private double ms2kts;
    private double NM2m;
    private double hp2W;
    private double g0;
    private double Per2c;
    private double missValue;
    private double rho_water;
    private double nu_water;
    private Integer int_precision;
    private double twentyfour;
    private double twelve;
    private double N2kN;

    public Const(){
        this.deg2rad = Math.PI/180;
        this.m2ft= 3.2808399;
        this.ms2kts=1.94384449;
        this.NM2m=1852.00;
        this.hp2W=0.7457*Math.pow(10,3);
        this.g0=9.80665;
        this.Per2c=(this.ms2kts*9.8)/(2*Math.PI);
        this.missValue=-99;
        this.rho_water=1029;
        this.nu_water=1.004*Math.pow(10,-6);
        this.int_precision = Integer.MAX_VALUE;//uint precision
        this.twentyfour=24;
        this.twelve=12;
        this.N2kN=Math.pow(10,-3);
    }


}
