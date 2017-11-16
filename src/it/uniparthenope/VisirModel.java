package it.uniparthenope;

public class VisirModel {
    //input data
    private int bar_flag;
    private int timedep_flag;
    private Const constants;
    private Optim optim;
    private Ship ship;
    private Visualization visualization;
    private SpatialGrid sGrid;
    private TemporalGrid tGrid;
    private ExtremePoints extreme_pts;
    private DepartureParameters dep_datetime;
    private SafetyParameters safety;

    //approssimating computation time
    private long startTime;
    private long estimatedTime;


    //Constructor
    public VisirModel(int bar_flag, int timedep_flag){//Initialize with standard values defined in settings.m
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

    public int getBar_flag() {
        return bar_flag;
    }

    public int getTimedep_flag() {
        return timedep_flag;
    }

    public long getEstimatedTime() {
        return estimatedTime;
    }

    /******************************************************/
}
