package it.uniparthenope;

public class SpatialGrid {
    private long numberOfNeighbors;
    private long numberOfShortestPaths;
    private double nodesLargeN;
    private long invStepFields;
    private double inv_step;
    private long minNoGridPoints;
    private double DB_bbox__lat_max;
    private double DB_bbox__lon_min;
    private double DB_bbox__lat_min;
    private double DB_bbox__lon_max;
    private double bbox__lat_max;
    private double bbox__lon_min;
    private double bbox__lat_min;
    private double bbox__lon_max;

    public SpatialGrid(){
        this.numberOfNeighbors = 24;
        this.numberOfShortestPaths = 1;
        this.nodesLargeN = 1.0 * Math.pow(10,5);//Limited by ram availability
        this.invStepFields = 4;//inverse grid size [deg^{-1}] of original wind forecast fields
        this.minNoGridPoints = 4; //minimum number of grid points along either direction (meridional or zonal) - TWA inconsistencies may result from tto small such a value
    }

    public void setInvStepFields(long flag){//relocatable model input fields
        if(flag == 2){
            this.invStepFields = 60; //inverse grid size [deg^{-1}]
        }
    }

    public void setInv_step(double inv_step){
        this.inv_step = inv_step;
    }

    public void setDB_bbox__lat_max(double DB_bbox__lat_max) {
        this.DB_bbox__lat_max = DB_bbox__lat_max;
    }

    public void setDB_bbox__lon_min(double DB_bbox__lon_min) {
        this.DB_bbox__lon_min = DB_bbox__lon_min;
    }

    public void setDB_bbox__lat_min(double DB_bbox__lat_min) {
        this.DB_bbox__lat_min = DB_bbox__lat_min;
    }

    public void setDB_bbox__lon_max(double DB_bbox__lon_max) {
        this.DB_bbox__lon_max = DB_bbox__lon_max;
    }

    public void setBbox__lat_max(double bbox__lat_max) {
        this.bbox__lat_max = bbox__lat_max;
    }

    public void setBbox__lon_min(double bbox__lon_min) {
        this.bbox__lon_min = bbox__lon_min;
    }

    public void setBbox__lat_min(double bbox__lat_min) {
        this.bbox__lat_min = bbox__lat_min;
    }

    public void setBbox__lon_max(double bbox__lon_max) {
        this.bbox__lon_max = bbox__lon_max;
    }

    public double getNodesLargeN() {
        return nodesLargeN;
    }

    public long getInvStepFields() {
        return invStepFields;
    }

    public long getMinNoGridPoints() {
        return minNoGridPoints;
    }

    public long getNumberOfNeighbors() {
        return numberOfNeighbors;
    }

    public long getNumberOfShortestPaths() {
        return numberOfShortestPaths;
    }

    public double getDB_bbox__lat_max() {
        return DB_bbox__lat_max;
    }

    public double getDB_bbox__lon_min() {
        return DB_bbox__lon_min;
    }

    public double getDB_bbox__lat_min() {
        return DB_bbox__lat_min;
    }

    public double getDB_bbox__lon_max() {
        return DB_bbox__lon_max;
    }

    public double getBbox__lat_max() {
        return bbox__lat_max;
    }

    public double getBbox__lon_min() {
        return bbox__lon_min;
    }

    public double getBbox__lat_min() {
        return bbox__lat_min;
    }

    public double getBbox__lon_max() {
        return bbox__lon_max;
    }

    public double getInv_step(){
        return  this.inv_step;
    }
}
