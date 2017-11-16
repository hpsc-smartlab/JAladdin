package it.uniparthenope;

public class Optim {
    //from settings.m
    private long relocAnlsCode;//code of the "analysis" of the  relocatable wave model

    private long intentional_speed_red;
    private long windModel;
    private long waveModel;

    public Optim(){
        this.relocAnlsCode = 20;
    }

    public void OptimizationParameters(){
        //getting data from optim_pars.json parsing
        MyJSONParser parser = new MyJSONParser("optim_pars.json");
        this.intentional_speed_red = parser.getValueAsLong("intentional_speed_red");
        this.windModel = parser.getValueAsLong("windModel");
        this.waveModel = parser.getValueAsLong("waveModel");
    }

    public long getRelocAnlsCode() {
        return relocAnlsCode;
    }

    public long getIntentional_speed_red() {
        return intentional_speed_red;
    }

    public long getWindModel() {
        return windModel;
    }

    public long getWaveModel() {
        return waveModel;
    }

}
