package it.uniparthenope;

import it.uniparthenope.Debug.MyFileWriter;
import it.uniparthenope.Parser.MyCSVParser;
import it.uniparthenope.Parser.MyJSONParser;

import java.io.FileInputStream;
import java.io.InputStream;
import java.io.Serializable;
import java.util.ArrayList;

public class ExtremePoints implements Serializable {
    //Route extreme points
    private double start_lat;
    private double start_lon;
    private double end_lat;
    private double end_lon;
    private double bbox_deltaLat_u;
    private double bbox_deltaLon_l;
    private double bbox_deltaLat_d;
    private double bbox_deltaLon_r;
    private long minCoastDist;
    private double UL_lat;
    private double UL_lon;
    private double DR_lat;
    private double DR_lon;
    private String cycType;
    private double pseudoG;


    public ExtremePoints(){
        //getting data from extrema_pars.json parsing
        MyJSONParser parser = new MyJSONParser("extrema_pars.json");
        this.start_lat = parser.getValueAsDouble("start_lat");
        this.start_lon = parser.getValueAsDouble("start_lon");
        this.end_lat = parser.getValueAsDouble("end_lat");
        this.end_lon = parser.getValueAsDouble("end_lon");
        this.bbox_deltaLat_u = parser.getValueAsDouble("bbox_deltaLat_u");
        this.bbox_deltaLon_l = parser.getValueAsDouble("bbox_deltaLon_l");
        this.bbox_deltaLat_d = parser.getValueAsDouble("bbox_deltaLat_d");
        this.bbox_deltaLon_r = parser.getValueAsDouble("bbox_deltaLon_r");
        this.minCoastDist = parser.getValueAsLong("minCoastDist");
        if(!RegionCheck()){
            System.out.println("Extreme points integrity check violated.");
            MyFileWriter debug = new MyFileWriter("","debug",false);
            debug.WriteLog("ExtremePoints constructor: Extreme points integrity check violated.");
            debug.CloseFile();
            System.exit(0);
        }
        this.cycType = "dd"; //cycloid type: 'id' % 'inverted_descent'; 'dd' % 'direct_descent'
    }

    private boolean RegionCheck(){
        // check if start or end points are in or out of permitted region
        boolean integrityCheck = true;
        try{
            InputStream stream = new FileInputStream("inputFiles/graph/freeedges_DB.dat.csv");
            MyCSVParser parser = new MyCSVParser(stream);
            ArrayList<Double> values = parser.getDoubleValues();
            this.UL_lat = values.get(0);
            this.UL_lon = values.get(1);
            this.DR_lat = values.get(2);
            this.DR_lon = values.get(3);
            if( ( this.start_lat > (this.UL_lat - this.bbox_deltaLat_u) ) ||
                    ( this.start_lon < (this.UL_lon - this.bbox_deltaLon_l) ) ||
                    ( this.end_lat < (this.DR_lat - this.bbox_deltaLat_d) ) ||
                    ( this.end_lon > (this.DR_lon - this.bbox_deltaLon_r) ) ){
                integrityCheck = false;
            }
        } catch(Exception e){
            integrityCheck = false;
            e.printStackTrace();
        }
        return integrityCheck;
    }

    public double getStart_lat() {
        return start_lat;
    }

    public double getStart_lon() {
        return start_lon;
    }

    public double getEnd_lat() {
        return end_lat;
    }

    public double getEnd_lon() {
        return end_lon;
    }

    public double getBbox_deltaLat_u() {
        return bbox_deltaLat_u;
    }

    public double getBbox_deltaLon_l() {
        return bbox_deltaLon_l;
    }

    public double getBbox_deltaLat_d() {
        return bbox_deltaLat_d;
    }

    public double getBbox_deltaLon_r() {
        return bbox_deltaLon_r;
    }

    public long getMinCoastDist() {
        return minCoastDist;
    }

    public double getUL_lat() {
        return UL_lat;
    }

    public double getUL_lon() {
        return UL_lon;
    }

    public double getDR_lat() {
        return DR_lat;
    }

    public double getDR_lon() {
        return DR_lon;
    }

    public void setStart_lat(double start_lat) {
        this.start_lat = start_lat;
    }

    public void setStart_lon(double start_lon) {
        this.start_lon = start_lon;
    }

    public void setEnd_lat(double end_lat) {
        this.end_lat = end_lat;
    }

    public void setEnd_lon(double end_lon) {
        this.end_lon = end_lon;
    }

    public void setCycType(String cycType) {
        this.cycType = cycType;
    }

    public void setBbox_deltaLat_u(double bbox_deltaLat_u) {
        this.bbox_deltaLat_u = bbox_deltaLat_u;
    }

    public void setBbox_deltaLon_l(double bbox_deltaLon_l) {
        this.bbox_deltaLon_l = bbox_deltaLon_l;
    }

    public void setBbox_deltaLat_d(double bbox_deltaLat_d) {
        this.bbox_deltaLat_d = bbox_deltaLat_d;
    }

    public void setBbox_deltaLon_r(double bbox_deltaLon_r) {
        this.bbox_deltaLon_r = bbox_deltaLon_r;
    }

    public String getCycType() {
        return cycType;
    }

    public double getPseudoG() {
        return pseudoG;
    }

    public void setPseudoG(double pseudoG) {
        this.pseudoG = pseudoG;
    }
}
