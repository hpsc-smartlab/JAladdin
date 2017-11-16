package it.uniparthenope;

public class ExtremePoints {
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

}
