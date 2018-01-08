package it.uniparthenope;

import it.uniparthenope.Parser.MyJSONParser;

import java.io.Serializable;

public class Visualization implements Serializable {
    private long resamp_factor;
    private long env_fields;
    private long waypoint_info;
    private long safegram;
    private long shipSpeedLUT;
    private char encounter_TF;
    private long scientific_mode;
    private long edgedelay_timedep;
    private long FnTilde_fitPlot;
    private long waveDisp;
    private long bathy;
    private long epsOutput;
    private long H_s;
    private long T_w;
    private long lambda;
    private long wind;
    private long graphData;
    private double customMax;
    private long fig_x;
    private long env_forcing;


    public Visualization(){
        this.resamp_factor = 24;//from settings.m. GDM paper (can be = 8. in this case, reference to IEEE paper)
    }

    public void VisualizationParameters(){
        //getting data from visualization_pars.json parsing
        MyJSONParser parser = new MyJSONParser("visualization_pars.json");
        this.env_fields = parser.getValueAsLong("env_fields");
        this.waypoint_info = parser.getValueAsLong("waypoint_info");
        this.safegram = parser.getValueAsLong("safegram");
    }

    public void setData(){
        this.shipSpeedLUT = 1;
        this.encounter_TF = 'T'; //Encounter period ('T') or frequency ('F')
        this.scientific_mode = 0; //0: bathy lsm after sea-over-land; % 1: no bathy lsm
        this.edgedelay_timedep = 0;
        this.FnTilde_fitPlot = 0;
        this.waveDisp = 0;
        this.bathy = 0;
        this.epsOutput = 1;
        this.H_s = 1; //Significant wave height % GMD, route 1924_i (GMD case study #1)
        this.T_w = 0; //Peak wave period % GMD, route 2067_e (GMD case study #3)
        this.lambda = 0; //Wavelength   % GMD, route 1589_c (GMD case study #2)
        this.wind = 1; //wind
        this.graphData = 0; //save graph to disk
        this.customMax = 3.5; //forcing max value for colorbar scale - throught different fields (wHeight,wPeriod,wLength)
        this.fig_x = 900;
        this.env_forcing = 0;
    }

    public long getResamp_factor() {
        return resamp_factor;
    }

    public long getEnv_fields() {
        return env_fields;
    }

    public long getWaypoint_info() {
        return waypoint_info;
    }

    public long getSafegram() {
        return safegram;
    }

    public long getShipSpeedLUT() {
        return shipSpeedLUT;
    }

    public char getEncounter_TF() {
        return encounter_TF;
    }

    public long getScientific_mode() {
        return scientific_mode;
    }

    public long getEdgedelay_timedep() {
        return edgedelay_timedep;
    }

    public long getFnTilde_fitPlot() {
        return FnTilde_fitPlot;
    }

    public long getWaveDisp() {
        return waveDisp;
    }

    public long getBathy() {
        return bathy;
    }

    public long getEpsOutput() {
        return epsOutput;
    }

    public long getH_s() {
        return H_s;
    }

    public long getT_w() {
        return T_w;
    }

    public long getLambda() {
        return lambda;
    }

    public long getWind() {
        return wind;
    }

    public long getGraphData() {
        return graphData;
    }

    public double getCustomMax() {
        return customMax;
    }

    public long getFig_x() {
        return fig_x;
    }

    public long getEnv_forcing(){
        return env_forcing;
    }

    public void setWaypoint_info(long waypoint_info) {
        this.waypoint_info = waypoint_info;
    }

    public void setSafegram(long safegram) {
        this.safegram = safegram;
    }

    public void setScientific_mode(long scientific_mode) {
        this.scientific_mode = scientific_mode;
    }

    public void setH_s(long h_s) {
        H_s = h_s;
    }

    public void setCustomMax(double customMax) {
        this.customMax = customMax;
    }

    public void setEnv_forcing(long env_forcing) {
        this.env_forcing = env_forcing;
    }
}
