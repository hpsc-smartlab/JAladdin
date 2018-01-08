package it.uniparthenope;


import java.io.FileInputStream;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

//NOTE:
//visualization (in settings.m file) class not provided in this version.
public class Main {

    public static void main(String[] args) {
        VisirModel v = new VisirModel(2,2, 1);//Same input provided for testing .m files
        v.Start();
    }
}
