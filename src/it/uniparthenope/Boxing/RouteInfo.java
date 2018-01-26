package it.uniparthenope.Boxing;

public class RouteInfo {
    Dijkstra2DResults route;
    double computationTime;//in seconds;

    public RouteInfo(Dijkstra2DResults route, double computationTime) {
        this.route = route;
        this.computationTime = computationTime;
    }

    public Dijkstra2DResults getRoute() {
        return route;
    }

    public double getComputationTime() {
        return computationTime;
    }
}
