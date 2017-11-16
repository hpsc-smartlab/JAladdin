package it.uniparthenope;

public class SafetyParameters {
    private long par_rolling;
    private long pureLossStab;
    private long surfRiding;

    public SafetyParameters(){
        //getting data from safety_pars.json parsing
        MyJSONParser parser = new MyJSONParser("safety_pars.json");
        this.par_rolling = parser.getValueAsLong("par_rolling");
        this.pureLossStab = parser.getValueAsLong("pureLossStab");
        this.surfRiding = parser.getValueAsLong("surfRiding");
    }

    public long getPar_rolling() {
        return par_rolling;
    }

    public long getPureLossStab() {
        return pureLossStab;
    }

    public long getSurfRiding() {
        return surfRiding;
    }
}
