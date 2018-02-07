package it.uniparthenope.Boxing;

import java.util.LinkedList;

public class RouteInfo {
    private int type;
    private double computationTime;//in seconds;
    private double[] partialTimes;
    private double[] Dr_cum;
    private double[] spd;
    private double[] dt_out;
    private double[] sh_course;
    private double[] theta_VMC;
    private double cost;
    private LinkedList<Integer> path;



    public RouteInfo(Dijkstra2DResults route, double computationTime, double[] partialTimes, get_Info_at_NodeResults routeInfo, int type) {
        //this.route = route;
        this.cost = route.getCost();
        this.path = route.getPath();
        this.computationTime = computationTime;
        this.partialTimes = partialTimes;
        Dr_cum = routeInfo.getDr_cum();
        spd = routeInfo.getV_opt();
        dt_out = routeInfo.getEdge_delay();
        sh_course = routeInfo.getTheta_opt();
        theta_VMC = routeInfo.getTheta_VMC();
        this.type = type;
    }

    public RouteInfo(DijkstraTimeResults route, double computationTime, get_Info_at_NodeResults routeInfo, int type){
        this.cost = route.getCost();
        this.path = route.getPath();
        this.partialTimes = route.getPartial_times();
        this.type = type;
        this.computationTime = computationTime;
        Dr_cum = routeInfo.getDr_cum();
        spd = routeInfo.getV_opt();
        dt_out = routeInfo.getEdge_delay();
        sh_course = routeInfo.getTheta_opt();
        theta_VMC = routeInfo.getTheta_VMC();
    }


    public double getComputationTime() {
        return computationTime;
    }

    public double[] getPartialTimes(){
        return partialTimes;
    }

    public double[] getDr_cum() {
        return Dr_cum;
    }

    public double[] getSpd() {
        return spd;
    }

    public double[] getDt_out() {
        return dt_out;
    }

    public double[] getSh_course() {
        return sh_course;
    }

    public double[] getTheta_VMC() {
        return theta_VMC;
    }

    public int getType() {
        return type;
    }

    public double getCost() {
        return cost;
    }

    public LinkedList<Integer> getPath() {
        return path;
    }

}
