package it.uniparthenope;


import java.io.FileInputStream;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

//NOTE:
//visualization (in settings.m file) class not provided in this version.
public class Main {

    public static void main(String[] args) {
        VisirModel v = new VisirModel(2,2);//Same input provided for testing .m files
        v.LoadData();
//        try {
//            InputStream stream = new FileInputStream("inputFiles/graph/freeedges_DB.dat.csv");
//            MyCSVParser parser = new MyCSVParser(stream);
//            ArrayList<Double> values = parser.getDoubleValues();
//            for(Double value: values){
//                System.out.println("Element: "+value);
//            }
//
//        } catch (Exception e){
//            e.printStackTrace();
//        }
        System.out.println("Data loaded");
    }
}
