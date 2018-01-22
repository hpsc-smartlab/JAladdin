package it.uniparthenope;

import java.util.ArrayList;

public class DangerIndexes implements Cloneable {
    private ArrayList<Integer> resonance1;
    private ArrayList<Integer> resonance2;
    private ArrayList<Integer> surfRiding;
    private ArrayList<Integer> pureLossStab;

    public DangerIndexes() {
        resonance1 = new ArrayList<>();
        resonance2 = new ArrayList<>();
        surfRiding = new ArrayList<>();
        pureLossStab = new ArrayList<>();
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        Object obj = super.clone();

        DangerIndexes idx = (DangerIndexes) obj;

        //deep cloning
        idx.setResonance1(resonance1);
        idx.setResonance2(resonance2);
        idx.setSurfRiding(surfRiding);
        idx.setPureLossStab(pureLossStab);
        return idx;
    }

    public ArrayList<Integer> getResonance1() {
        return resonance1;
    }

    public ArrayList<Integer> getResonance2() {
        return resonance2;
    }

    public void setResonance1(ArrayList<Integer> resonance1) {
        this.resonance1 = resonance1;
    }

    public void setResonance2(ArrayList<Integer> resonance2) {
        this.resonance2 = resonance2;
    }

    public ArrayList<Integer> getSurfRiding() {
        return surfRiding;
    }

    public void setSurfRiding(ArrayList<Integer> surfRiding) {
        this.surfRiding = surfRiding;
    }

    public ArrayList<Integer> getPureLossStab() {
        return pureLossStab;
    }

    public void setPureLossStab(ArrayList<Integer> pureLossStab) {
        this.pureLossStab = pureLossStab;
    }

}
