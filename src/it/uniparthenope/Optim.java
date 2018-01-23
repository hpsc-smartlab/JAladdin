package it.uniparthenope;

import it.uniparthenope.Parser.JSONManager;
import it.uniparthenope.Parser.MyJSONParser;
import org.json.simple.parser.ParseException;

import java.io.IOException;


public class Optim {
    //from settings.m
    private long relocAnlsCode;//code of the "analysis" of the  relocatable wave model

    private long intentional_speed_red;
    private long windModel;
    private long waveModel;

    public Optim(){
        this.relocAnlsCode = 20;
    }

    public Optim(boolean flag) throws IOException, ParseException {
        if(flag==true){
            JSONManager reader = new JSONManager();
            reader.initReading("SerializedObjects/Optim.json");
            relocAnlsCode = reader.retrieveLong("relocAnlsCode");
            intentional_speed_red = reader.retrieveLong("intentional_speed_red");
            windModel = reader.retrieveLong("windModel");
            waveModel = reader.retrieveLong("waveModel");
            reader.dispose();
        }
    }

    public void OptimizationParameters(){
        //getting data from optim_pars.json parsing
        MyJSONParser parser = new MyJSONParser("optim_pars.json");
        this.intentional_speed_red = parser.getValueAsLong("intentional_speed_red");
        this.windModel = parser.getValueAsLong("windModel");
        this.waveModel = parser.getValueAsLong("waveModel");
    }


    public void saveState() throws IOException {
        JSONManager writer = new JSONManager();
        writer.initWriting("SerializedObjects/Optim.json");
        writer.putLong("relocAnlsCode", relocAnlsCode);
        writer.putLong("intentional_speed_res", intentional_speed_red);
        writer.putLong("windModel", windModel);
        writer.putLong("waveModel", waveModel);
        writer.dispose();
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
