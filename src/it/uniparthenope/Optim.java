package it.uniparthenope;

public class Optim {
    //from settings.m
    private int relocAnlsCode;//code of the "analysis" of the  relocatable wave model

    private int intentional_speed_red;
    private int windModel;
    private int waveModel;

    public Optim(){
        this.relocAnlsCode = 20;
    }

    public void OptimizationParameters(){
        //getting data from optim_pars.json parsing
        MyJSONParser parser = new MyJSONParser("optim_pars.json");
        this.intentional_speed_red = parser.getValueAsInt("intentional_speed_red");
        this.windModel = parser.getValueAsInt("windModel");
        this.waveModel = parser.getValueAsInt("waveModel");
    }

    public int getRelocAnlsCode() {
        return relocAnlsCode;
    }

    public int getIntentional_speed_red() {
        return intentional_speed_red;
    }

    public int getWindModel() {
        return windModel;
    }

    public int getWaveModel() {
        return waveModel;
    }

}
