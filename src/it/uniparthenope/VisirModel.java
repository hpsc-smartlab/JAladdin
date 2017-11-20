package it.uniparthenope;

import java.util.ArrayList;
import java.util.Collections;

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
    private Double[][] vel_LUT;
    private ArrayList<Double> H_array_m;
    private Double minWind;
    private Double maxWind;
    private Double polar_twa;
    private Double polar_tws;

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

    private class Output{
        public ArrayList<Double> H_array_m;
        public Double[][] vel_LUT;
    }

    private void ship_Model(){//called from vessel_Response.m, ship_model.m implementation
        //Look-up table for involuntary ship speed reduction in waves.
        this.ship.setFn_max(this.constants.getMs2kts(), this.constants.getG0());
        //LUT independent variables:
        long Nh = 25; //40
        long Nl = 40; //30
        double Hmax = 8.0; //[m]
        long ih2 = (long) Math.floor(Nh/2);
        long il2 = (long) Math.floor(Nl/2);

        //Significant wave height
        //Adding to H_array_m (Nh-1) elements between 10^-1 and 10^log10(Hmax)
        ArrayList<Double> H_array_m = Utility.logspace(-1.0, Hmax, Nh-1);
        H_array_m.add(0, 0.0);//adding 0 at first element of H_array_m
        //preallocations:
        Double[][] vel_LUT = Utility.zeros((int) Nh+1, this.ship.getP_level_hp().size()+1);
        Double[][] Rc_LUT = Utility.zeros((int) Nh+1, this.ship.getP_level_hp().size()+1);
        Double[][] Raw_LUT = Utility.zeros((int) Nh+1, this.ship.getP_level_hp().size()+1);
        ArrayList<Double> P_level_thro = new ArrayList<Double>();
        Double max = Collections.max(this.ship.getP_level_hp());
        for(Double element : this.ship.getP_level_hp()){
            P_level_thro.add((100*element) / max);
        }
        //pars for v_Bowditch:
        Double m2ft = 3.2808399;
        Double a3_ref = 0.0248; //[kts/ ft^2]
        Double a2_ref = 0.0165; //[kts/ ft^2]
        Double a1_ref = 0.0083; //[kts/ ft^2]
        Double v0_ref = 18.0; //kts
        //LUT computation
        for(int ih = 1; ih< Nh; ih++){
            this.ship.ship_resistance(H_array_m.get(ih), this.constants);
            for(int j = 1; j<this.ship.getP_level_hp().size(); j++){
                vel_LUT[ih][j] = this.ship.getV_out_kts().get(j);
                Rc_LUT[ih][j] = this.ship.getR_c().get(j);
                Raw_LUT[ih][j] = this.ship.getR_aw().get(j);
                //v_Bowditch(ih)= max( 0, vel_LUT(1,1)/v0_ref* ( v0_ref - a2_ref *  (m2ft*H_array_m(ih)).^2) );
            }
        }
        System.out.println("5");
        this.ship.ship_resistance(0.0, this.constants);
        ArrayList<Double> v0 = this.ship.getV_out_kts();
        ArrayList<Double> Rc0 = this.ship.getR_c();
        ArrayList<Double> Raw0 = this.ship.getR_aw();
//        for ip=1: numel(P_level_hp)
//          disp([ min(vel_LUT(:,ip)), max(vel_LUT(:,ip)) ])
//        end
        //graphical output not implemented.
        this.H_array_m = H_array_m;
        this.vel_LUT = vel_LUT;
    }

    public void vessel_Response(){//vessel_Response.m implementation
        this.ship.shipPowerLevels();
        this.ship_Model();
        this.maxWind = Double.NaN;
        this.minWind = Double.NaN;
        this.polar_twa = Double.NaN;
        this.polar_tws = Double.NaN;
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

    public Double[][] getVel_LUT() {
        return vel_LUT;
    }

    public ArrayList<Double> getH_array_m() {
        return H_array_m;
    }

    /******************************************************/
}
