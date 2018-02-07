package it.uniparthenope.ShipModel;

import it.uniparthenope.Boxing.ship_ModelResults;
import it.uniparthenope.Const;
import it.uniparthenope.Debug.MyFileWriter;
import it.uniparthenope.Parser.MyJSONParser;
import it.uniparthenope.Utility;

import java.util.ArrayList;
import java.util.Collections;

public class Motorboat implements Ship {

    /*defined in settings.m*/
    private long sailType;
    private double etaEngine;
    /***********************/
    private long vessType;
    private long sailClass;
    private long P_max_hp;
    private double maxv;
    private double length;
    private double beam;
    private double draught;
    private double obs_roll_period;
    private double roll_period;
    private double bottomDraughtfactor;
    private double TEsaturFactor;
    private double Rspectrum;
    private double Nvel;
    private ArrayList<Double> P_level_hp;
    private double Fn_max;
    private ArrayList<Double> v_out_kts;
    private ArrayList<Double> R_c;
    private ArrayList<Double> R_aw;
    private double maxWind;
    private double minWind;
    private double min_lambda_L;
    private double max_lambda_L;
    private double[][] ship_v_LUT;
    private ArrayList<Double> H_array_m;
    private ArrayList<Double> polar_twa;
    private ArrayList<Double> polar_tws;

    public Motorboat(){
        this.sailType = 4;
        this.etaEngine = 0.7;
        this.ship_v_LUT = new double[0][0];
        this.H_array_m = new ArrayList<>();
        this.polar_twa = new ArrayList<>();
        this.polar_tws = new ArrayList<>();
    }

    public void LoadVesselParameters(String path, String outDir){
        //getting data from ship_pars.json parsing
        MyJSONParser parser = new MyJSONParser(path);
        this.vessType = parser.getValueAsLong("vessType");
        this.sailClass = parser.getValueAsLong("sailClass");
        this.P_max_hp = parser.getValueAsLong("P_max_hp");
        this.maxv = parser.getValueAsDouble("maxv");
        this.length = parser.getValueAsDouble("length(m)");
        this.beam = parser.getValueAsDouble("beam(m)");
        this.draught = parser.getValueAsDouble("draught(m)");
        this.obs_roll_period = parser.getValueAsDouble("obs_roll_period(s)");
        if(!IntegrityCheck()){
            System.out.println("Ship integrity check violated.");
            MyFileWriter debug = new MyFileWriter("","debug",false, outDir);
            debug.WriteLog("LoadVesselParameters: Ship integrity check violated.");
            debug.CloseFile();
            System.exit(0);
        }
        if( this.vessType != this.sailType){
            this.roll_period = this.obs_roll_period;
        } else {
            this.roll_period = 1;
        }
        this.bottomDraughtfactor = 1; //bottomDraughtfactor*draught is threshold for effect of sea bottom on enhanced frictional
        this.TEsaturFactor = 10; //max allowed TE in plots is ship.roll_period*ship.TEsaturFactor
        this.Rspectrum = 0;
        this.Nvel = 1; //Default value for steps in power reduction.
        this.P_level_hp = new ArrayList<Double>();
    }

    private boolean IntegrityCheck(){//namelist_check.m file porting
        boolean check = true;
        if(this.vessType != this.sailType){
            if( this.P_max_hp <= 0 ||
                    this.maxv <= 0 ||
                    this.length <= 0 ||
                    this.beam <=0 ||
                    this.draught <= 0 ||
                    this.obs_roll_period <= 0){
                check = false;
            }
        }
        return check;
    }

    public void setStepsInPowerReduction(long flag){
        if(flag==1){
            this.Nvel = 7;
        } else { this.Nvel = 1; }
        if (this.vessType == this.sailType){ //sailboat
            this.Nvel = 1;
        }
    }

    public void shipPowerLevels(){
        if(this.Nvel ==1){
            //constant power level:
            this.P_level_hp.add(1.0);
        } else{
            this.P_level_hp = Utility.linspace(0.1, 1.0, (long) this.Nvel);
            Collections.reverse(this.P_level_hp);
        }
        for(int i=0;i<this.P_level_hp.size();i++){
            this.P_level_hp.set(i,this.P_level_hp.get(i)*this.P_max_hp);
        }
    }

    @Override
    public void shipResistance(double wHeight, Const constants) {
        //ship speed computation out of wave height, wavelenght and towing Power.
        //results are nearly linear in power for not too large w.height

        //Ship parameters

        //http://www.tadroberts.ca/services/small-boats/displacement/fidler19
        //Double TT = 0.5; //m
        //Double BB = 2.1; //m
        //Double Lwl = 5.5; //m
        //Double Delta = 1.2 * Math.pow(10,3); //kg
        //Double v_max_kts = 5.0; //kts
        //this.P_max_hp = 40;

        //Galatea
        //http://www.marina.difesa.it/uominimezzi/navi/Pagine/Galatea.aspx > scheda tecnica
        //Double TT = 2.5; //m
        //Double BB = 12.6; //m
        //Double Lwl = 39.2; //m
        //Double Delta = 415.0 * Math.pow(10,3); //kg
        //Double v_max_kts = 13.0; //kts
        //this.P_max_hp = 1872;

        double Lwl = this.length;
        double BB = this.beam;
        double TT = this.draught;
        double Swet = 2*BB*TT + BB*Lwl + 2*TT*Lwl; //parallelepiped hull

        //Conversions:
        double twopi = 2* Math.PI;
        ArrayList<Double> P_w = new ArrayList<Double>();
        for (double element : this.P_level_hp){
            P_w.add(constants.getHp2W() * element);
        }
        double P_max_w = constants.getHp2W() * this.P_max_hp;
        double v_max_ms = this.maxv / constants.getMs2kts();
        //COMPUTATIONS:
        //wave added resistance: P=c2ca* v^2
        double xi = 0.5 * wHeight; //xi= wave amplitude = 1/2 wave height
//         Double sigma_aw_tilde = 1.0;
//         Double sigma_aw_tilde = 3.5; //VISIR fishing vessel
//         Double sigma_aw_tilde = 2.2; //VISIR ferry
//         Double sigma_aw_tilde = cos(cost.deg2rad*alpha/2) .* sigma_aw_tilde
        double sigma_aw_tilde = this.ship_Raw_alexandersson(Lwl,BB,TT);
        double Fn_tilde = this.ship_Fn_tilde(constants);
        double k2 = sigma_aw_tilde * constants.getRho_water() * Math.sqrt( constants.getG0()/Math.pow(Lwl,3) ) * Math.pow(xi,2) * Math.pow(BB,2) /
                ( 2* this.etaEngine * Fn_tilde);
        // Double alpha = 0.0; //angle of encounter in deg from North
        // k2ca = k2 * cos(constants.deg2rad * alpha);

        //Calm water resistance: P = k3 * v^3
        double k3 = P_max_w / Math.pow(v_max_ms,3); //W*s^3/m^3
        long Rr_exp = 2; // default is = 2; Referee #1 cp. also Harvald Eq.4.3.27
        long n_exp = Rr_exp - 2;
        if(n_exp < 0){
            System.out.println("Residual resistance decreasing with vessel speed!");
        }
        double v_search = v_max_ms;
        this.v_out_kts = new ArrayList<Double>();
        this.R_c = new ArrayList<Double>();
        this.R_aw = new ArrayList<Double>();
        this.v_out_kts.add(0,Double.NaN);
        this.R_c.add(0,Double.NaN);
        this.R_aw.add(0,Double.NaN);
        String num_method ="Newton";
        for(int ip=0; ip<P_w.size(); ip++){
            double k0 = -P_w.get(ip);
            double v_out_ms=0;
            if(num_method == "Newton") {
                v_out_ms = Utility.Newton(k3, k2, k0, n_exp, v_max_ms, v_search);//used as fzero matlab approssimation
            } else {
                v_out_ms = Utility.nr_cubic(k2/k3, 0 , k0/k3); //nr_cubic generic method
            }
            v_search = v_out_ms;
            this.v_out_kts.add(ip,constants.getMs2kts()*v_out_ms);

            //Resistances:
            double val = constants.getN2kN() * this.etaEngine * k3 * Math.pow(v_out_ms,2+n_exp) / Math.pow(v_max_ms, n_exp); //kN = P_c/v
            R_c.add(ip, val);
            double val2 = constants.getN2kN() * this.etaEngine * k2 * Math.pow(v_max_ms, 1); //kN = P_aw/v
            R_aw.add(ip, val2);
        }
    }
    private double ship_Raw_alexandersson(double L, double B, double T){//from ship_Raw_alexandersson.m file
        // Peak value of added wave resistance according to
        // nonlinear regression of model simulations
        // by Alexandersson 2009, Master Thesis, KTH Centre for Naval Architecture

        // % Alexandersson ship_1:
        // % L=142.5;
        // % B=22.5;
        // % T=8;
        // % Fn=.21;
        // % cp=.57;
        // % Raw_peak=10.8

        // % VISIR fishing vessel:
        // % L=16;
        // % B=5.2;
        // % T=1.6;
        // % Fn=.21;
        // % Raw_peak= 6.5

        // % VISIR ferry:
        // % L=86.6;
        // % B=24;
        // % T=3.4;
        // % Fn=.21;
        // % Raw_peak= 4.3

        // % % fishing vessel C.484 (INSEAN 1962)
        // % L=14.25;
        // % B=4.52;
        // % T=1.9;
        // % cp=.6;
        // % Fn=.21;
        // % % Raw_peak= 8

        double Bn = B/L;
        double Tn = T/L;

        //Original formulation by Alexandersson 2009:
        // % Raw_peak= 206.6 *  L_CG^0.6376 * Bn^(-1.2121) * Tn^ 0.6247 * kyy^1.3611 * Fn^aexp;
        // % I assume that L_CG is a mistake and that L_CGn is meant instead!
        //
        // % My assumptions:
        // % L_CGn=.50;  % = Lcg/L  (Lcg = Longitudinal center gravity)
        // % kyy=.25;    % = Ryy/L  (Ryy= Radius of pitch gyration)
        // % Cp^0.0440    = 1;
        // % L_CGn^0.6376 = 0.642781
        // % kyy^1.3611   = 0.1515;
        //
        // % Resulting expression:
        double Raw_peak = 20.0 * Math.pow(Bn,-1.21) * Math.pow(Tn,0.62);
        return Raw_peak;
    }

    private double ship_Fn_tilde(Const constant){//from ship_Fn_tilde.m
        //identification of linear slope fitting the Alexandersson dependence on Fn in the region
        //of values defined by [Fn_max/speedRedFact, Fn_max]
        // L_test is ship length in m
        // v_test is ship speed in kts
        // g0=9.80665; % m/s^2
        // ms2kts=1.94384449;
        //
        // lop
        // v_test =  22 ; % kts
        // L_test =  86 ; % m
        double L_test = this.length;
        double v_test = this.maxv;
        double Fn_max = v_test / (constant.getMs2kts()*Math.sqrt(constant.getG0()*L_test));
        long Ndots = 100;
        long speedRedFact = 1000;
        ArrayList<Double> fn = Utility.linspace(Fn_max/speedRedFact, Fn_max, Ndots);
        ArrayList<Double> alex = new ArrayList<Double>();
        for(double element : fn){
            alex.add(Math.pow(element,0.6377)); //Alexandersson
        }
        double slope = linfit_origin(fn,alex);
        double fn_tilde = 1/slope;
        return fn_tilde;

    }

    private double linfit_origin(ArrayList<Double> x, ArrayList<Double> y){
        // Least square fit of linear funztion through origin.
        // cp. Numerical Recipes,
        // http://hebb.mit.edu/courses/9.29/2002/readings/c15-2.pdf
        // set a=0 in Eq.(15.2.6), with sigma=1 in definitions Eq.(15.2.4)
        // and post-processing after setting a=0
        int S = x.size();
        double Sx = 0.0;
        ArrayList<Double> Sx2 = new ArrayList<Double>();
        for(double element : x){
            Sx += element;
            Sx2.add(Math.pow(element,2));
        }
        double Sy = 0.0;
        for(Double element : y){
            Sy += element;
        }
        double Sxx = 0.0;
        for(Double element : Sx2){
            Sxx += element;
        }
        double Delta = S*Sxx - Math.pow(Sx,2);
        return (S*Sxx*Sy/Sx - Sx*Sy)/Delta;
    }

    @Override
    public void vesselResponse(Const constants, String outdir) {
        MyFileWriter logFile =new MyFileWriter("","",true, outdir);
        logFile.WriteLog("Ship model response (motorboat)...");
        logFile.CloseFile();
        this.shipPowerLevels();
        ship_Model(constants);
        this.maxWind = Double.NaN;
        this.minWind = Double.NaN;
        this.polar_tws = new ArrayList<>();
        this.polar_twa = new ArrayList<>();
    }

    private void ship_Model(Const constants){//called from vessel_Response.m, ship_model.m implementation
        //Look-up table for involuntary ship speed reduction in waves.
        this.setFn_max(constants.getMs2kts(), constants.getG0());
        //LUT independent variables:
        long Nh = 25; //40
        long Nl = 40; //30
        double Hmax = 8.0; //[m]
        //long ih2 = (long) Math.floor(Nh/2);
        //long il2 = (long) Math.floor(Nl/2);

        //Significant wave height
        //Adding to H_array_m (Nh-1) elements between 10^-1 and 10^log10(Hmax)
        ArrayList<Double> H_array_m = Utility.logspace(-1.0, Math.log10(Hmax), Nh-1);
        H_array_m.add(0, 0.0);//adding 0 at first element of H_array_m

        //preallocations:
        double[][] vel_LUT = new double[(int)Nh][this.getP_level_hp().size()];

        double max = Collections.max(this.getP_level_hp());

        //LUT computation
        for(int ih = 0; ih<Nh; ++ih){
            this.shipResistance(H_array_m.get(ih), constants);
            for(int j = 0; j<this.getP_level_hp().size(); ++j){
                vel_LUT[ih][j] = this.getV_out_kts().get(j);
            }
        }
        this.shipResistance(0.0, constants);

        this.ship_v_LUT = vel_LUT;
        this.H_array_m = H_array_m;
    }



    /***GETTER AND SETTER***/

    public long getSailType() {
        return sailType;
    }

    public void setSailType(long sailType) {
        this.sailType = sailType;
    }

    public double getEtaEngine() {
        return etaEngine;
    }

    public void setEtaEngine(double etaEngine) {
        this.etaEngine = etaEngine;
    }

    public long getVessType() {
        return vessType;
    }

    public void setVessType(long vessType) {
        this.vessType = vessType;
    }

    public long getSailClass() {
        return sailClass;
    }

    public void setSailClass(long sailClass) {
        this.sailClass = sailClass;
    }

    public long getP_max_hp() {
        return P_max_hp;
    }

    public void setP_max_hp(long p_max_hp) {
        P_max_hp = p_max_hp;
    }

    public double getMaxv() {
        return maxv;
    }

    public void setMaxv(double maxv) {
        this.maxv = maxv;
    }

    public double getLength() {
        return length;
    }

    public void setLength(double length) {
        this.length = length;
    }

    public double getBeam() {
        return beam;
    }

    public void setBeam(double beam) {
        this.beam = beam;
    }

    public double getDraught() {
        return draught;
    }

    public void setDraught(double draught) {
        this.draught = draught;
    }

    public double getObs_roll_period() {
        return obs_roll_period;
    }

    public void setObs_roll_period(double obs_roll_period) {
        this.obs_roll_period = obs_roll_period;
    }

    public double getRoll_period() {
        return roll_period;
    }

    public void setRoll_period(double roll_period) {
        this.roll_period = roll_period;
    }

    public double getBottomDraughtfactor() {
        return bottomDraughtfactor;
    }

    public void setBottomDraughtfactor(double bottomDraughtfactor) {
        this.bottomDraughtfactor = bottomDraughtfactor;
    }

    public double getTEsaturFactor() {
        return TEsaturFactor;
    }

    public void setTEsaturFactor(double TEsaturFactor) {
        this.TEsaturFactor = TEsaturFactor;
    }

    public double getRspectrum() {
        return Rspectrum;
    }

    public void setRspectrum(double rspectrum) {
        Rspectrum = rspectrum;
    }

    public double getNvel() {
        return Nvel;
    }

    public void setNvel(double nvel) {
        Nvel = nvel;
    }

    public ArrayList<Double> getP_level_hp() {
        return P_level_hp;
    }

    public void setP_level_hp(ArrayList<Double> p_level_hp) {
        P_level_hp = p_level_hp;
    }

    public double getFn_max() {
        return Fn_max;
    }

    public void setFn_max(double ms2kts, double g0){
        this.Fn_max = this.maxv/(ms2kts * Math.sqrt(g0*this.length));
    }

    public ArrayList<Double> getV_out_kts() {
        return v_out_kts;
    }

    public void setV_out_kts(ArrayList<Double> v_out_kts) {
        this.v_out_kts = v_out_kts;
    }

    public ArrayList<Double> getR_c() {
        return R_c;
    }

    public void setR_c(ArrayList<Double> r_c) {
        R_c = r_c;
    }

    public ArrayList<Double> getR_aw() {
        return R_aw;
    }

    public void setR_aw(ArrayList<Double> r_aw) {
        R_aw = r_aw;
    }

    public double getMaxWind() {
        return maxWind;
    }

    public void setMaxWind(double maxWind) {
        this.maxWind = maxWind;
    }

    public double getMinWind() {
        return minWind;
    }

    public void setMinWind(double minWind) {
        this.minWind = minWind;
    }

    public double getMin_lambda_L() {
        return min_lambda_L;
    }

    public void setMin_lambda_L(double min_lambda_L) {
        this.min_lambda_L = min_lambda_L;
    }

    public double getMax_lambda_L() {
        return max_lambda_L;
    }

    public void setMax_lambda_L(double max_lambda_L) {
        this.max_lambda_L = max_lambda_L;
    }

    public double[][] getShip_v_LUT() {
        return ship_v_LUT;
    }

    public void setShip_v_LUT(double[][] ship_v_LUT) {
        this.ship_v_LUT = ship_v_LUT;
    }

    public ArrayList<Double> getH_array_m() {
        return H_array_m;
    }

    public void setH_array_m(ArrayList<Double> h_array_m) {
        H_array_m = h_array_m;
    }

    public ArrayList<Double> getPolar_twa() {
        return polar_twa;
    }

    public void setPolar_twa(ArrayList<Double> polar_twa) {
        this.polar_twa = polar_twa;
    }

    public ArrayList<Double> getPolar_tws() {
        return polar_tws;
    }

    public void setPolar_tws(ArrayList<Double> polar_tws) {
        this.polar_tws = polar_tws;
    }
}
