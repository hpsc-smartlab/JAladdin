package it.uniparthenope;


import java.io.FileInputStream;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

//NOTE:
//visualization (in settings.m file) class not provided in this version.
public class Main {

    public static void main(String[] args) {
        //mode=0: normal execution without serialize objects
        //mode=1: normal execution serializing objects
        //mode=2: debugging mode: load data deserializing objects
        VisirModel v = new VisirModel(2,2, 2);//Same input provided for testing .m files
        v.Start();
    }
}
