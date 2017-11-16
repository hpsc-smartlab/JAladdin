package it.uniparthenope;

//NOTE:
//visualization (in settings.m file) class not provided in this version.
public class Main {

    public static void main(String[] args) {
        VisirModel v = new VisirModel(2,2);//Same input provided for testing .m files
        v.LoadData();
        System.out.println("Data loaded");
    }
}
