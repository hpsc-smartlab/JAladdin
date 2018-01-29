package it.uniparthenope.Boxing;

public class RouteInfo {
    private Dijkstra2DResults route;
    private int type;
    private double computationTime;//in seconds;
    private double[] partialTimes;
    private double[] Dr_cum;
    private double[] spd;
    private double[] dt_out;
    private double[] sh_course;
    private double[] theta_VMC;


    public RouteInfo(Dijkstra2DResults route, double computationTime, double[] partialTimes, get_Info_at_NodeResults routeInfo, int type) {
        this.route = route;
        this.computationTime = computationTime;
        this.partialTimes = partialTimes;
        Dr_cum = routeInfo.getDr_cum();
        spd = routeInfo.getV_opt();
        dt_out = routeInfo.getEdge_delay();
        sh_course = routeInfo.getTheta_opt();
        theta_VMC = routeInfo.getTheta_VMC();
        this.type = type;
    }

    public Dijkstra2DResults getRoute() {
        return route;
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
}
