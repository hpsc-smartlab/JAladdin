package it.uniparthenope;

public class SafetyParameters {
    private int par_rolling;
    private int pureLossStab;
    private int surfRiding;

    public SafetyParameters(){
        //getting data from safety_pars.json parsing
        MyJSONParser parser = new MyJSONParser("safety_pars.json");
        this.par_rolling = parser.getValueAsInt("par_rolling");
        this.pureLossStab = parser.getValueAsInt("pureLossStab");
        this.surfRiding = parser.getValueAsInt("surfRiding");
    }

    public int getPar_rolling() {
        return par_rolling;
    }

    public int getPureLossStab() {
        return pureLossStab;
    }

    public int getSurfRiding() {
        return surfRiding;
    }
}
