package it.uniparthenope;

public class Visualization {
    private int resamp_factor;
    private int env_fields;
    private int waypoint_info;
    private int safegram;


    public Visualization(){
        this.resamp_factor = 24;//from settings.m. GDM paper (can be = 8. in this case, reference to IEEE paper)
    }

    public void VisualizationParameters(){
        //getting data from visualization_pars.json parsing
        MyJSONParser parser = new MyJSONParser("visualization_pars.json");
        this.env_fields = parser.getValueAsInt("env_fields");
        this.waypoint_info = parser.getValueAsInt("waypoint_info");
        this.safegram = parser.getValueAsInt("safegram");
    }

    public int getResamp_factor() {
        return resamp_factor;
    }

    public int getEnv_fields() {
        return env_fields;
    }

    public int getWaypoint_info() {
        return waypoint_info;
    }

    public int getSafegram() {
        return safegram;
    }
}
