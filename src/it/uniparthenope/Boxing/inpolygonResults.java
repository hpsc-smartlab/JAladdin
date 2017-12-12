package it.uniparthenope.Boxing;

public class inpolygonResults {
    private boolean in;
    private boolean on;

    public inpolygonResults(){
        this.in = false;
        this.on = false;
    }

    public inpolygonResults(boolean in, boolean on){
        this.in = in;
        this.on = on;
    }

    public boolean getIn() {
        return in;
    }

    public void setIn(boolean in) {
        this.in = in;
    }

    public boolean getOn() {
        return on;
    }

    public void setOn(boolean on) {
        this.on = on;
    }
}
