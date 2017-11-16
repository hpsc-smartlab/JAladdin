package it.uniparthenope;

public class Visualization {
    private long resamp_factor;
    private long env_fields;
    private long waypoint_info;
    private long safegram;


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
}
