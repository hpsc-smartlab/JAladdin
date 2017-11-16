package it.uniparthenope;

public class VisirModel {
    //input data
    private long bar_flag;
    private long timedep_flag;
    private Const constants;
    private Optim optim;
    private Ship ship;
    private Visualization visualization;
    private SpatialGrid sGrid;
    private TemporalGrid tGrid;
    private ExtremePoints extreme_pts;
    private DepartureParameters dep_datetime;
    private SafetyParameters safety;
    private EnvironmentalFields forcing;

    //approssimating computation time
    private long startTime;
    private long estimatedTime;


    //Constructor
    public VisirModel(long bar_flag, long timedep_flag){//Initialize with standard values defined in settings.m
        this.Tic();//Start the "timer"

        //bar_flag = 1: fresh compution of edges not crossing coastline (mode 1 of GMD-D paper)
        //bar_flag = 2: edges not crossing coastline read out from DB   (mode 2 of GMD-D paper)
        if((bar_flag!=1) && (bar_flag!=2)) { //If the input incorrect, set bar_flag = 2 as default value.
            this.bar_flag = 2;
        } else {
            this.bar_flag = bar_flag;
        }
        // timedep_flag=0: static algorithm with fields at time step #1
        // timedep_flag=2: time-dependent method
        if((timedep_flag!=0) && (timedep_flag!=2)) { //Same thing for timedep_flag.
            this.timedep_flag = 2;
        } else {
            this.timedep_flag = timedep_flag;
        }
        this.constants = new Const();
        this.optim = new Optim();
        this.ship = new Ship();
        this.sGrid = new SpatialGrid();
        this.tGrid = new TemporalGrid();
        this.visualization = new Visualization();

    }

    public void LoadData(){//Loading data parsing them from json file defined in inputFiles folder
        this.extreme_pts = new ExtremePoints();
        this.dep_datetime = new DepartureParameters();
        this.ship.LoadVesselParameters();
        this.safety = new SafetyParameters();
        this.optim.OptimizationParameters();
        this.visualization.VisualizationParameters();
    }

    public void CalculateParameters(){//Set some parameters based on other parameters
        this.ship.setStepsInPowerReduction(this.optim.getIntentional_speed_red());
        this.sGrid.setInvStepFields(this.optim.getWaveModel());
        this.forcing = new EnvironmentalFields(this.ship.getVessType(), this.ship.getSailType(),this.optim.getWindModel());
        this.safety.setCriteria();
        this.visualization.setData();
        if(this.forcing.getAnalytic() == 1){
            this.extreme_pts.setStart_lat(41.0);
            this.extreme_pts.setStart_lon(10.0);
            this.extreme_pts.setEnd_lat(40.5);
            this.extreme_pts.setEnd_lon(11.0);
            this.extreme_pts.setCycType("id");
            this.dep_datetime.setYear(14);
            this.dep_datetime.setMonth(1);
            this.dep_datetime.setDay(1);
            this.dep_datetime.setHour(0);
            this.dep_datetime.setMin(0);
            this.visualization.setEnv_forcing(1);
            this.visualization.setWaypoint_info(0);
            this.visualization.setSafegram(0);
            this.visualization.setScientific_mode(1);
            this.visualization.setH_s(1);
            this.visualization.setCustomMax(Double.NaN);
            this.ship.setVessType(1);
            this.timedep_flag = 0;
        }
    }


    private void Tic(){//Equivalent to MATLAB tic function
        this.startTime = System.nanoTime();
    }

    private void Toc(){
        this.estimatedTime = System.nanoTime() - this.startTime;
    }

    /*****************Getter methods*********************/
    public Const getConstants() {
        return this.constants;
    }

    public Optim getOptim() {
        return this.optim;
    }

    public Ship getShip() {
        return this.ship;
    }

    public SpatialGrid getsGrid() {
        return this.sGrid;
    }

    public TemporalGrid gettGrid() {
        return this.tGrid;
    }

    public ExtremePoints getExtreme_pts() {
        return extreme_pts;
    }

    public SafetyParameters getSafety() {
        return safety;
    }

    public DepartureParameters getDep_datetime() {
        return dep_datetime;
    }

    public Visualization getVisualization() {
        return visualization;
    }

    public long getBar_flag() {
        return bar_flag;
    }

    public long getTimedep_flag() {
        return timedep_flag;
    }

    public long getEstimatedTime() {
        return estimatedTime;
    }

    /******************************************************/
}
