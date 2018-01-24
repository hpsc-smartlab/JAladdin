package it.uniparthenope;

import it.uniparthenope.Parser.JSONManager;
import org.json.simple.parser.ParseException;

import java.io.IOException;

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
    private double DB_xi;
    private double DB_yi;
    private double DB_xf;
    private double DB_yf;
    private long DB_Nx;
    private long DB_Ny;
    private double xi;
    private double yi;
    private double xf;
    private double yf;
    private double inset_area;
    private long inset_Nx;
    private long inset_Ny;
    private long freenodes;
    private double min_start_dist;
    private double min_end_dist;
    private long node_start;
    private long node_end;
    private double node_start_lat;
    private double node_start_lon;
    private double node_end_lat;
    private double node_end_lon;
    private double theta_gdt;


    public SpatialGrid(){
        this.numberOfNeighbors = 24;
        this.numberOfShortestPaths = 1;
        this.nodesLargeN = 1.0 * Math.pow(10,5);//Limited by ram availability
        this.invStepFields = 4;//inverse grid size [deg^{-1}] of original wind forecast fields
        this.minNoGridPoints = 4; //minimum number of grid points along either direction (meridional or zonal) - TWA inconsistencies may result from tto small such a value
    }

    public SpatialGrid(boolean flag) throws IOException, ParseException, java.text.ParseException{
        if(flag==true){
            JSONManager reader = new JSONManager();
            reader.initReading("SerializedObjects/SpatialGrid.json");
            numberOfNeighbors = reader.retrieveLong("numberOfNeighbors");
            numberOfShortestPaths = reader.retrieveLong("numberOfShortestPaths");
            nodesLargeN = reader.retrieveDouble("nodesLargeN");
            invStepFields = reader.retrieveLong("invStepFields");
            inv_step = reader.retrieveDouble("inv_step");
            minNoGridPoints = reader.retrieveLong("minNoGridPoints");
            DB_bbox__lat_max = reader.retrieveDouble("DB_bbox__lat_max");
            DB_bbox__lat_min = reader.retrieveDouble("DB_bbox__lat_min");
            DB_bbox__lon_max = reader.retrieveDouble("DB_bbox__lon_max");
            DB_bbox__lon_min = reader.retrieveDouble("DB_bbox__lon_min");
            bbox__lat_max = reader.retrieveDouble("bbox__lat_max");
            bbox__lat_min = reader.retrieveDouble("bbox__lat_min");
            bbox__lon_max = reader.retrieveDouble("bbox__lon_max");
            bbox__lon_min = reader.retrieveDouble("bbox__lon_min");
            DB_xi = reader.retrieveDouble("DB_xi");
            DB_yi = reader.retrieveDouble("DB_yi");
            DB_xf = reader.retrieveDouble("DB_xf");
            DB_yf = reader.retrieveDouble("DB_yf");
            DB_Nx = reader.retrieveLong("DB_Nx");
            DB_Ny = reader.retrieveLong("DB_Ny");
            xi = reader.retrieveDouble("xi");
            yi = reader.retrieveDouble("yi");
            xf = reader.retrieveDouble("xf");
            yf = reader.retrieveDouble("yf");
            inset_area = reader.retrieveDouble("inset_area");
            inset_Nx = reader.retrieveLong("inset_Nx");
            inset_Ny = reader.retrieveLong("inset_Ny");
            freenodes = reader.retrieveLong("freenodes");
            min_start_dist = reader.retrieveDouble("min_start_dist");
            min_end_dist = reader.retrieveDouble("min_end_dist");
            node_start = reader.retrieveLong("node_start");
            node_end = reader.retrieveLong("node_end");
            node_start_lat = reader.retrieveDouble("node_start_lat");
            node_end_lat = reader.retrieveDouble("node_end_lat");
            node_start_lon = reader.retrieveDouble("node_start_lon");
            node_end_lon = reader.retrieveDouble("node_end_lon");
            theta_gdt = reader.retrieveDouble("theta_gdt");
            reader.dispose();
        }
    }

    public void saveState() throws IOException, java.text.ParseException{
        JSONManager writer = new JSONManager();
        writer.initWriting("SerializedObjects/SpatialGrid.json");
        writer.putLong("numberOfNeighbors", numberOfNeighbors);
        writer.putLong("numberOfShortestPaths",numberOfShortestPaths);
        writer.putDouble("nodesLargeN", nodesLargeN);
        writer.putLong("invStepFields", invStepFields);
        writer.putDouble("inv_step", inv_step);
        writer.putLong("minNoGridPoints", minNoGridPoints);
        writer.putDouble("DB_bbox__lat_max", DB_bbox__lat_max);
        writer.putDouble("DB_bbox__lat_min", DB_bbox__lat_min);
        writer.putDouble("DB_bbox__lon_max", DB_bbox__lon_max);
        writer.putDouble("DB_bbox__lon_min", DB_bbox__lon_min);
        writer.putDouble("bbox__lat_max", bbox__lat_max);
        writer.putDouble("bbox__lat_min", bbox__lat_min);
        writer.putDouble("bbox__lon_max", bbox__lon_max);
        writer.putDouble("bbox__lon_min", bbox__lon_min);
        writer.putDouble("DB_xi", DB_xi);
        writer.putDouble("DB_yi", DB_yi);
        writer.putDouble("DB_xf", DB_xf);
        writer.putDouble("DB_yf", DB_yf);
        writer.putLong("DB_Nx", DB_Nx);
        writer.putLong("DB_Ny", DB_Ny);
        writer.putDouble("xi", xi);
        writer.putDouble("yi", yi);
        writer.putDouble("xf", xf);
        writer.putDouble("yf", yf);
        writer.putDouble("inset_area", inset_area);
        writer.putLong("inset_Nx", inset_Nx);
        writer.putLong("inset_Ny", inset_Ny);
        writer.putLong("freenodes", freenodes);
        writer.putDouble("min_start_dist", min_start_dist);
        writer.putDouble("min_end_dist", min_end_dist);
        writer.putLong("node_start", node_start);
        writer.putLong("node_end", node_end);
        writer.putDouble("node_start_lon", node_start_lon);
        writer.putDouble("node_start_lat", node_start_lat);
        writer.putDouble("node_end_lon", node_end_lon);
        writer.putDouble("node_end_lat", node_end_lat);
        writer.putDouble("theta_gdt", theta_gdt);
        writer.dispose();
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

    public double getDB_xi() {
        return DB_xi;
    }

    public void setDB_xi(double DB_xi) {
        this.DB_xi = DB_xi;
    }

    public double getDB_yi() {
        return DB_yi;
    }

    public void setDB_yi(double DB_yi) {
        this.DB_yi = DB_yi;
    }

    public double getDB_xf() {
        return DB_xf;
    }

    public void setDB_xf(double DB_xf) {
        this.DB_xf = DB_xf;
    }

    public double getDB_yf() {
        return DB_yf;
    }

    public void setDB_yf(double DB_yf) {
        this.DB_yf = DB_yf;
    }

    public long getDB_Nx() {
        return DB_Nx;
    }

    public void setDB_Nx(long DB_Nx) {
        this.DB_Nx = DB_Nx;
    }

    public long getDB_Ny() {
        return DB_Ny;
    }

    public void setDB_Ny(long DB_Ny) {
        this.DB_Ny = DB_Ny;
    }

    public double getXi() {
        return xi;
    }

    public void setXi(double xi) {
        this.xi = xi;
    }

    public double getYi() {
        return yi;
    }

    public void setYi(double yi) {
        this.yi = yi;
    }

    public double getXf() {
        return xf;
    }

    public void setXf(double xf) {
        this.xf = xf;
    }

    public double getYf() {
        return yf;
    }

    public void setYf(double yf) {
        this.yf = yf;
    }

    public double getInset_area() {
        return inset_area;
    }

    public void setInset_area(double inset_area) {
        this.inset_area = inset_area;
    }

    public long getInset_Nx() {
        return inset_Nx;
    }

    public void setInset_Nx(long inset_Nx) {
        this.inset_Nx = inset_Nx;
    }

    public long getInset_Ny() {
        return inset_Ny;
    }

    public void setInset_Ny(long inset_Ny) {
        this.inset_Ny = inset_Ny;
    }

    public long getFreenodes() {
        return freenodes;
    }

    public void setFreenodes(long freenodes) {
        this.freenodes = freenodes;
    }

    public double getMin_start_dist() {
        return min_start_dist;
    }

    public void setMin_start_dist(double min_start_dist) {
        this.min_start_dist = min_start_dist;
    }

    public double getMin_end_dist() {
        return min_end_dist;
    }

    public void setMin_end_dist(double min_end_dist) {
        this.min_end_dist = min_end_dist;
    }

    public long getNode_start() {
        return node_start;
    }

    public void setNode_start(long node_start) {
        this.node_start = node_start;
    }

    public long getNode_end() {
        return node_end;
    }

    public void setNode_end(long node_end) {
        this.node_end = node_end;
    }

    public double getNode_start_lat() {
        return node_start_lat;
    }

    public void setNode_start_lat(double node_start_lat) {
        this.node_start_lat = node_start_lat;
    }

    public double getNode_start_lon() {
        return node_start_lon;
    }

    public void setNode_start_lon(double node_start_lon) {
        this.node_start_lon = node_start_lon;
    }

    public double getTheta_gdt() {
        return theta_gdt;
    }

    public void setTheta_gdt(double theta_gdt) {
        this.theta_gdt = theta_gdt;
    }

    public double getNode_end_lat() {
        return node_end_lat;
    }

    public void setNode_end_lat(double node_end_lat) {
        this.node_end_lat = node_end_lat;
    }

    public double getNode_end_lon() {
        return node_end_lon;
    }

    public void setNode_end_lon(double node_end_lon) {
        this.node_end_lon = node_end_lon;
    }
}
