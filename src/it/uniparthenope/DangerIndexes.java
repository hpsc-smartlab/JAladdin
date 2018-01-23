package it.uniparthenope;

import it.uniparthenope.Parser.JSONManager;
import org.json.simple.parser.ParseException;

import java.io.IOException;
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

    public DangerIndexes(int idx) throws IOException, ParseException {
        JSONManager reader = new JSONManager();
        reader.initReading("SerializedObjects/DangerIndexes"+idx+".json");
        resonance1 = reader.retrieveIntegerArrayList("resonance1");
        resonance2 = reader.retrieveIntegerArrayList("resonance2");
        surfRiding = reader.retrieveIntegerArrayList("surfRiding");
        pureLossStab = reader.retrieveIntegerArrayList("pureLossStab");
        reader.dispose();
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

    public void saveState(int idx) throws IOException {
        JSONManager writer = new JSONManager();
        writer.initWriting("SerializedObjects/DangerIndexes"+idx+".json");
        writer.putIntegerArrayList("resonance1", resonance1);
        writer.putIntegerArrayList("resonance2", resonance2);
        writer.putIntegerArrayList("surfRiding", surfRiding);
        writer.putIntegerArrayList("pureLossStab", pureLossStab);
        writer.dispose();
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
