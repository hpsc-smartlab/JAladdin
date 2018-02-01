package it.uniparthenope;

public class InputPaths {
    private String dep_parameters;
    private String extr_parameters;
    private String optim_parameters;
    private String safety_parameters;
    private String ship_parameters;
    private String visualization_parameters;
    private String outDir;
    private String freeNodesDB;
    private String coastlineDB;
    private String bathymetryDB;
    private String analysisDB;
    private String forecastFile;
    private String freeEdgesDB;

    public InputPaths(){
        //default files
        this.outDir = "Output/";//default output directory
        this.dep_parameters = "inputFiles/datetime_pars.json";
        this.extr_parameters = "inputFiles/extrema_pars.json";
        this.optim_parameters = "inputFiles/optim_pars.json";
        this.safety_parameters = "inputFiles/safety_pars.json";
        this.ship_parameters = "inputFiles/ship_pars.json";
        this.visualization_parameters = "inputFiles/visualization_pars.json";
        this.freeNodesDB = "inputFiles/graph/freeNodes_DB.dat";
        this.coastlineDB = "inputFiles/coast/medf.map";
        this.bathymetryDB = "inputFiles/bathy/MedOneMin/med_one_min_single_vs2.nc";
        this.analysisDB = "inputFiles/fields/an_dates_DB.txt";
//        this.forecastFile = "inputFiles/wave/WW3/forecast/start__20150329.nc";
//        this.forecastFile = "inputFiles/wave/WW3/forecast/20180128.nc";
        this.forecastFile = "inputFiles/wave/WW3/forecast/start__20131226.nc";
        this.freeEdgesDB = "inputFiles/graph/freeedges_DB.dat";
    }

    public InputPaths(String dep_parameters,String extr_parameters,String optim_parameters,
                      String safety_parameters, String ship_parameters, String visualization_parameters,
                      String outDir, String freeNodesDB, String coastlineDB, String bathymetryDB,
                      String analysisDB, String forecastFile, String freeEdgesDB){
        if(outDir==null)
            this.outDir = "Output/";//default output directory
        else
            this.outDir = outDir+"/";
        if(dep_parameters==null)
            this.dep_parameters = "inputFiles/datetime_pars.json";
        else
            this.dep_parameters = dep_parameters;
        if(extr_parameters==null)
            this.extr_parameters = "inputFiles/extrema_pars.json";
        else
            this.extr_parameters = extr_parameters;
        if(optim_parameters==null)
            this.optim_parameters = "inputFiles/optim_pars.json";
        else
            this.optim_parameters = optim_parameters;
        if(safety_parameters==null)
            this.safety_parameters = "inputFiles/safety_pars.json";
        else
            this.safety_parameters = safety_parameters;
        if(ship_parameters==null)
            this.ship_parameters = "inputFiles/ship_pars.json";
        else
            this.ship_parameters = ship_parameters;
        if(visualization_parameters==null)
            this.visualization_parameters = "inputFiles/visualization_pars.json";
        else
            this.visualization_parameters = visualization_parameters;
        if(freeNodesDB==null)
            this.freeNodesDB = "inputFiles/graph/freeNodes_DB.dat";
        else
            this.freeNodesDB = freeNodesDB;
        if(coastlineDB==null)
            this.coastlineDB = "inputFiles/coast/medf.map";
        else
            this.coastlineDB = coastlineDB;
        if(bathymetryDB==null)
            this.bathymetryDB = "inputFiles/bathy/MedOneMin/med_one_min_single_vs2.nc";
        else
            this.bathymetryDB = bathymetryDB;
        if(analysisDB==null)
            this.analysisDB = "inputFiles/fields/an_dates_DB.txt";
        else
            this.analysisDB = analysisDB;
        if(forecastFile==null)
            this.forecastFile = "inputFiles/wave/WW3/forecast/start__20150329.nc";
        else
            this.forecastFile = forecastFile;
        if(freeEdgesDB==null)
            this.freeEdgesDB = "inputFiles/graph/freeedges_DB.dat";
        else
            this.freeEdgesDB = freeEdgesDB;

    }

    public String getDep_parameters() {
        return dep_parameters;
    }

    public String getExtr_parameters() {
        return extr_parameters;
    }

    public String getOptim_parameters() {
        return optim_parameters;
    }

    public String getSafety_parameters() {
        return safety_parameters;
    }

    public String getShip_parameters() {
        return ship_parameters;
    }

    public String getVisualization_parameters() {
        return visualization_parameters;
    }

    public String getOutDir() {
        if(outDir==null)
            return "Output/";
        else
            return outDir;
    }

    public String getFreeNodesDB() {
        return freeNodesDB;
    }

    public String getCoastlineDB() {
        return coastlineDB;
    }

    public String getBathymetryDB() {
        return bathymetryDB;
    }

    public String getAnalysisDB() {
        return analysisDB;
    }

    public String getForecastFile() {
        return forecastFile;
    }

    public String getFreeEdgesDB() {
        return freeEdgesDB;
    }
}
