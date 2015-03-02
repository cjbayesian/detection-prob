#include<iostream>
#include<fstream>
#include<cstdlib>
#include<math.h>
#include<cmath>
#include<string>
#include<sstream>
#include<modRank.h>
#include<vbc_vector.h>
#include<Simplex.h>
#include<metro_hastingsMD.h>
#include<iomanip>
#include<sys/stat.h>

#define MATHLIB_STANDALONE 1
#include <Rmath.h>

#ifndef GRAV_H
#define GRAV_H

using namespace std;
using namespace modRank;
using namespace vbc_lib;
using namespace simplex;
using namespace mcmcMD;

// GLOBALS //
    bool impose_management = FALSE;
    int man_dim = 1;
    int alpha_order = 1; //polynomial order of the env-alpha function //
    bool traf_mat_precompute = FALSE;
    bool show_ = FALSE;
    _vbc_vec<float> glb_true_params;
    float gen_pdet;
    bool include_undetected_sites=TRUE;
    bool compensate_force=FALSE;
    bool only_known=FALSE;
    bool force_invasions=FALSE;
    int glb_r;
    int rngseed=1234;
    bool one_time_traf=TRUE;
    float gen_prop=0;
    int run_type=2; //1:MLE, 2:MCMC, 3:sim spread from posterior
    int post_length=0;
    bool env=FALSE; //set to FALSE for reproducing Gertzen.
    bool fit_pdet=FALSE;
    int n_sim_for_smooth=100; //number of sims for smoothing the likelihood surface
    bool ll=FALSE;
    int est_env=13;
    int est_pdet=1;
    bool sim=FALSE;
    bool hetero_pdet=FALSE;
    bool gridded=TRUE;
    bool boot=TRUE;
    float min_d = 0.01;
    float max_d = 5;
    float min_e = 0.01;
    float max_e = 3;
    int tmat_res = 5; //20
    int glb_traf_array_ind;
    float pdet=1; //detection probability to be fit
    float d_par=1; //0.288344;   //to be fit
    float e_par=1;        //to be fit
    float c_par=1;          //to be fit
    float gamma_par=0;      //to be fit
    float glb_alpha=0.00001;          //to be fit (as a function of chem_vars)
    int n_lakes=1646;
    int n_sources=213;      //4496; with non grided Oi
    int n_val_lakes;    //Should be 100 -- need to review data parsing code.
    int n_chem_var=13;
    int n_pdet_var=1;
    int from_year=1989;
    int to_year=2009;
    int n_sampled=0;
    float fixed_d=0;
    bool new_tr_par=TRUE;
    string parseeds_file("no_env.par");
    _vbc_vec<float> d_matrix;
    _vbc_vec<float> chem_pars(1,n_chem_var+1); //to be fit
    _vbc_vec<float> pdet_pars(1,n_pdet_var+1);
    _vbc_vec<float> traf_mat;
    _vbc_vec<float> traf_mat_array;
    _vbc_vec<int> t_vec;
    _vbc_vec<int> sampled_index;
    _vbc_vec<int> val_lakes_index(1,n_val_lakes);
    _vbc_vec<int> val_lakes_inv_index(1,n_lakes);
    _vbc_vec<float> theta_vec(1,n_lakes);
    _vbc_vec<float> phi_vec(1,n_lakes);
    _vbc_vec<string> parnames;
    string sites_file_name("../2010_bytho_data/lakes_processed_normalized.csv");
    string env_file_name("");
    string pod_pred_file_name("");
    string disc_file_name("");
    string dist_file_name("../2010_bytho_data/distance_matrix_grd.csv"); //full: "../2010_bytho_data/cj_roadconn.tab"
    string Oi_file_name("../2010_bytho_data/Oi_grd.csv"); //full: "../2010_bytho_data/Oi.csv"
    string output_folder("");
//debugging output //
    ofstream track_lake("output/tl.dat");
    ofstream t_file("output/t_mcmc.dat");
    ofstream l_file("output/l_mcmc.dat");
    ofstream par_file("output/par_mcmc.dat");
    ofstream alpha_mcmc_file("output/alpha_mcmc.dat");
    ofstream pr_inv_file("output/sim_prop_inv.csv");
// CLASSES //
class clsLake
	{
	public:
		float Wj,x,y,alpha,Uj;
        bool invaded;
        bool val_invaded;
        int discovered;
        int last_abs;
        int status_2010;
        int id;
        _vbc_vec<float> chem;
	    _vbc_vec<float> Uij;
	    _vbc_vec<float> Qjt;
        _vbc_vec<float> pp;
        _vbc_vec<float> PInv;
        _vbc_vec<float> PInv_inst;
        _vbc_vec<float> pdet_var;
	};
class clsSource
	{
	public:
		float Ai;
		int Oi;
        _vbc_vec<float> Xit; //prop_invaded
        _vbc_vec<float> Dij;
		_vbc_vec<float> Gij;
	};

class clsState
	{
	public:
        int year;
        int n_inv;
        int n_new_inv;
        int n_u_inv;
        int n_calc_pp;
        _vbc_vec<int> inv;
        _vbc_vec<int> new_inv;
		_vbc_vec<int> u_inv;
		_vbc_vec<int> calc_pp_index;    
        
	};

_vbc_vec<clsLake> lakes(1,n_lakes);
_vbc_vec<clsSource> sources(1,n_sources);
_vbc_vec<clsState> state(from_year,to_year);


void args(int ,char *);
void sample_t_l(int);
void which_sampled_or_valid();
void read_data();
void read_po_data();

float calc_pr_I(int);
float calc_alpha(int);
float calc_pdet(int);
void calc_pp();
void calc_pp_validation(_vbc_vec<int>);
void update_sim_pp(int,_vbc_vec<int>);
void update_pp_l_hood(int, int);
int which_min(int,int);
int which_max(int,int);
void write_t();
void write_par();
void write_traf_mat();
void calc_traf_mat_part(_vbc_vec<int>);
void clear_traf_mat();
void calc_traf_mat_array();
void write_inv_stat();
void write_pp();
float l_hood();
float l_hood_detp();
float l_hood_detp2();
float l_hood_detp3();
float l_hood_detp4();
float l_hood_detp5();
void mean_field();

float average(_vbc_vec<float> );
int wc_l(string);
int wc_column(string);
float add_log(float, float);
void normalize_vars();
void print_lake(int);

float pr_inv_not_det(int, int);

float MLE_l_hood(_vbc_vec<float> *,_vbc_vec<float> *);

void likelihood_wrapperMCMC(_vbc_vec<float> *, float *,int );
void likelihood_wrapperMCMC_MD(_vbc_vec<float> *, float *,int );
void likelihood_wrapperMCMC_MD_joint(_vbc_vec<float> *, float *,int );
float prior(float,int);
bool restrict_MCMC(float,int);

float prior_MD(_vbc_vec<float>,int);
bool restrict_MCMC_MD(_vbc_vec<float>);
_vbc_vec<float> predict_p(_vbc_vec<float>,_vbc_vec<int>,int); //pars,indicies
_vbc_vec<int> sample_w_replace(_vbc_vec<int>);
_vbc_vec<int> sample_wo_replace(_vbc_vec<int>,int);

void link_pars(_vbc_vec<string> ,_vbc_vec<float> );
void parse_par_seeds(string,_vbc_vec<float> *);

void save_beta_dists(_vbc_vec<float> * );
_vbc_vec<float> fit_beta(_vbc_vec<float> *,int, int);
float beta_lhood(_vbc_vec<float> * , _vbc_vec<float> *);
#endif;

////////////////////////////////////////////////////////////////////
void args(int argc,char *argv[])
{
	// decode arguments //
	int optind=1;
	while ((optind < argc) && (argv[optind][0]=='-')) 
	{
      string sw = argv[optind];
      if (sw=="-seeds") //Parameter seeds file (also indicates which params to fit)
		{    
			optind++;
            parseeds_file = argv[optind];
	   }
      else if (sw=="-runtype") //general flag for different run modes
      {
         optind++;
         run_type = atoi(argv[optind]);
      }
      else if (sw=="-rngseed") //rng seed integer
      {
         optind++;
         rngseed = atoi(argv[optind]);
      }
      else if (sw=="-sites") // path to file with sites data (id, location, area)
      {
         optind++;
         sites_file_name = argv[optind];
      }
      else if (sw=="-Oi") // path to file with source locations (dispersers)
      {
         optind++;
         Oi_file_name = argv[optind];
      }
      else if (sw=="-Env") // path to file with environmental covariariates
      {
         optind++;
         env_file_name = argv[optind];
      }
      else if (sw=="-PoDvars") // path to file with probability of detection predictors
      {
         optind++;
         pod_pred_file_name = argv[optind];
      }
      else if (sw=="-discs") // path to file with discoveries data
      {
         optind++;
         disc_file_name = argv[optind];
      }
      else if (sw=="-dist") // path to file with distances between sources and sites
      {
         optind++;
         dist_file_name = argv[optind];
      }
      else if (sw=="-of") // folder for outputs
      {
         optind++;
         output_folder = argv[optind];
      }
      else if (sw=="-ll") //-ll flag to call l_hood() with given pars
      {
            ll=TRUE; 
            sim = true;  
            optind++;         
            n_sources = atoi(argv[optind]);
            optind++;
            n_lakes = atoi(argv[optind]);
            optind++;
            d_par = atof(argv[optind]);
            optind++;
            e_par = atof(argv[optind]);
        }
        else if (sw=="-lld") //-ll flag to call l_hood() with given pars
        {    
            ll=TRUE; 
            sim = false;  
            optind++;         
            d_par = atof(argv[optind]);
            optind++;
            e_par = atof(argv[optind]);
        }
        else
		{
         cerr << "Unknown switch: " << argv[optind] << endl; 
		   	//<< "Options:\n\t-seeds <seed params file>  \n\t-runtype: <integer run type>\n" 
            //<< "Usage: ./grav -s <n_sources> <n_lakes>\n";
            exit(1);
		}
        optind++;
	} // end of arg-decoding
}

void inits()
{
   string tmp;
   ifstream init_file("inits.ini");

   init_file >> tmp;
   init_file >> post_length;
   init_file >> tmp;
   init_file >> n_sim_for_smooth;
   init_file >> tmp;
   init_file >> gridded;
   init_file >> tmp;
   init_file >> boot;
   init_file >> tmp;
   init_file >> fixed_d;
   init_file >> tmp;
   init_file >> to_year;
   init_file >> tmp;
   init_file >> tmat_res;

   init_file.close();
   if(fixed_d != 0)
      d_par = fixed_d;
}


void read_data()
{
    // Lakes and chemistry //
    ifstream l_file;
    l_file.open(sites_file_name.c_str()); 
    n_lakes = wc_l(sites_file_name.c_str());
    /* File structure:
    "Hectares",
    "UTM_X",
    "UTM_Y",
    "invaded",
    "B_DISC_YEAR",
    "LAST_ABS",
    "status_2010",
    "NAUT",
    "KKUT",
    "MGUT",
    "CAUT",
    "PPUT1",
    "SIO3UR",
    "DOC",
    "COLTR",
    "ALKTI",
    "ALKT",
    "PH",
    "COND25",
    "SECCHI.DEPTH"
    */
    traf_mat.redim(1,n_lakes,1,n_lakes);
    t_vec.redim(1,n_lakes);
    theta_vec.redim(1,n_lakes);
    phi_vec.redim(1,n_lakes);
    for(int i = 1;i<=n_lakes;i++)
	{
        lakes(i).PInv.redim(from_year,to_year);
        lakes(i).PInv_inst.redim(from_year+1,to_year);
        lakes(i).chem.redim(1,n_chem_var);
        lakes(i).pp.redim(from_year,to_year);
        lakes(i).pp(from_year)=0; //in li_hood() lakes invaded in 1989 don't contribute to likelihood.
        lakes(i).pdet_var.redim(1,n_pdet_var);

        l_file >> lakes(i).Wj;
        l_file >> lakes(i).x;
        l_file >> lakes(i).y;
        l_file >> lakes(i).invaded;
        l_file >> lakes(i).discovered;
        l_file >> lakes(i).last_abs;
        l_file >> lakes(i).status_2010;

        //Management vectors set to 1 (no cleaning stations)//
        theta_vec(i) = 1;
        phi_vec(i) = 1;

        // Catch for to_year < 2010.
        // modify lakes(i) to reflect missing information.
        /*
         lakes(i).val_invaded=0;
         if(lakes(i).discovered > 2009)
         {
            lakes(i).invaded=0;
            lakes(i).val_invaded=1;
            lakes(i).discovered=0;
            n_validation++;
            tmp_val_lakes_index(n_validation)=i;
            n_invaded_in_validation_set++;
            tmp_val_inv_lakes_index(n_invaded_in_validation_set)=i;
         }
         if(lakes(i).last_abs > 2009)
         {
            lakes(i).last_abs = 0;
            n_validation++;
            tmp_val_lakes_index(n_validation)=i;
         }

         */
        for(int j = 1;j<=n_chem_var;j++)
        {
            l_file >> lakes(i).chem(j);
        }
        for(int j = 1;j<=n_pdet_var;j++)
        {
            lakes(i).pdet_var(j) = runif(-2,2); // TODO Read from file if given //
        }

    }

    l_file.close();

    // Distance matrix //
    ifstream d_file;
    d_file.open(dist_file_name.c_str()); //distance_matrix.csv
    n_sources = wc_l(dist_file_name.c_str());
    
    d_matrix.redim(1,n_sources,1,n_lakes);
    sources.redim(1,n_sources);

    for(int i = 1;i<=n_sources;i++)
	 {
        sources(i).Dij.redim(1,n_lakes);
        sources(i).Gij.redim(1,n_lakes);
        sources(i).Xit.redim(from_year-1,to_year);
        sources(i).Xit(from_year-1)=0;
        for(int j = 1;j<=n_lakes;j++)
        {
            d_file >> sources(i).Dij(j);
        }
    }
    d_file.close();

    // Oi //
    ifstream o_file;
    o_file.open(Oi_file_name.c_str());

    for(int i = 1;i<=n_sources;i++)
	 {
         o_file >> sources(i).Oi;
         sources(i).Oi=sources(i).Oi;   
         //cout <<   sources(i).Oi << "\n";
    }
    o_file.close();

    for(int ch=1;ch<=n_chem_var+1;ch++)
      chem_pars(ch)=0; //chem_pars(ch)=runif(-0.01,0.01);

}
void read_po_data() //read in presence-only data
{

    /////////////////////////////////////////////////    
    // sites (lakes) //
    ifstream l_file;
    l_file.open(sites_file_name.c_str()); 
    n_lakes = wc_l(sites_file_name.c_str());
    traf_mat.redim(1,n_lakes,1,n_lakes);    
    t_vec.redim(1,n_lakes);
    theta_vec.redim(1,n_lakes);
    phi_vec.redim(1,n_lakes);

    if(wc_column(sites_file_name.c_str()) != 4){cerr << "**Sites file should contain 4 cols. <id> <x> <y> <area>\n";exit(1);}
    for(int i = 1;i<=n_lakes;i++)
	{
        lakes(i).chem.redim(1,n_chem_var);
        lakes(i).pdet_var.redim(1,n_pdet_var);
        lakes(i).discovered = 0;

        //Management vectors set to 1 (no cleaning stations)//
        theta_vec(i) = 1;
        phi_vec(i) = 1;

        l_file >> lakes(i).id;
        l_file >> lakes(i).x;
        l_file >> lakes(i).y;
        l_file >> lakes(i).Wj;

        //cerr << lakes(i).id << " " << lakes(i).x << " " << lakes(i).y << " " << lakes(i).Wj << "\n";
    }
    l_file.close();
    /////////////////////////////////////////////////


    /////////////////////////////////////////////////
    // Envirnomental Variabes //
    if(env_file_name == ""){cerr << "**Use the -Env flag to provide environmental covariariates file.\n";exit(1);}
    int tmp_n_lakes = wc_l(env_file_name.c_str());
    n_chem_var = wc_column(env_file_name.c_str()) - 1;
    cerr << "Reading " << n_chem_var << " env vars...\n";
    ifstream env_file;
    env_file.open(env_file_name.c_str());
    if(tmp_n_lakes != n_lakes){cerr << "**Environmental variables file not same length as Sites file**\n";exit(1);}
    for(int i = 1;i<=n_lakes;i++)
	{   
        lakes(i).chem.redim(1,n_chem_var);
        int tmp_id;
        env_file >> tmp_id;
        if(tmp_id != lakes(i).id){cerr << "Mismatch in environmental vars file. IDs not lined up with Sites file.\n";exit(1);}
        for(int j = 1;j<=n_chem_var;j++)
        {
            env_file >> lakes(i).chem(j);
        }
    }
    /////////////////////////////////////////////////


    /////////////////////////////////////////////////
    // Probability of detection predictors //
    if(pod_pred_file_name == ""){
        n_pdet_var = 0;
    } else {
        tmp_n_lakes = wc_l(pod_pred_file_name.c_str());
        n_pdet_var = wc_column(pod_pred_file_name.c_str());
        ifstream pod_pred_file;
        pod_pred_file.open(pod_pred_file_name.c_str());
        if(tmp_n_lakes != n_lakes){cerr << "**Environmental variables file not same length as Sites file**\n";}
        for(int i = 1;i<=n_lakes;i++)
	    {   
            lakes(i).pdet_var.redim(1,n_pdet_var);
            for(int j = 1;j<=n_pdet_var;j++)
            {
                pod_pred_file >> lakes(i).pdet_var(j); // TODO Read from file if given //
            }
        }
    }
    /////////////////////////////////////////////////
    

    /////////////////////////////////////////////////
    // Normalize predictors to zero mean, 1 std.
    normalize_vars();
    /////////////////////////////////////////////////



    /////////////////////////////////////////////////
    // Read in discoveries //
    if(disc_file_name == ""){
        cerr << "No Discovery data provided, simulating...\n";
    } else {
        from_year = 3000; //set in the future then update with earliest discovery //
        ifstream disc_file;
        disc_file.open(disc_file_name.c_str()); 
        int n_disc = wc_l(disc_file_name.c_str());
        if(wc_column(disc_file_name.c_str()) != 2){cerr << "Discoveries files should be two columns: <id> <year>\n";}
        if(n_disc == 0){cerr << "No discoveries provided\n";}
        int tmp_disc_id;
        int tmp_disc_year;
        for(int d = 1;d<=n_disc;d++)
	    {
            disc_file >> tmp_disc_id;
            disc_file >> tmp_disc_year;
            for(int i=1;i<=n_lakes;i++)
            {
                if(lakes(i).id == tmp_disc_id)
                {
                    lakes(i).discovered = tmp_disc_year;
                    cerr << lakes(i).id << " discovered in " << tmp_disc_year << "\n";
                    break;
                }
            }
            if(tmp_disc_year < from_year)
                from_year = tmp_disc_year;
        }
        disc_file.close();
    }
    /////////////////////////////////////////////////

    /////////////////////////////////////////////////
    // Init pp and mean field prob vectors //
    for(int i=1;i<=n_lakes;i++)
    {
        //if(lakes(i).id == 1668056) //LAKE ONTARIO (SEED)
        //    lakes(i).discovered == from_year;
        lakes(i).PInv.redim(from_year,to_year);
        lakes(i).PInv_inst.redim(from_year+1,to_year);
        lakes(i).pp.redim(from_year,to_year);
        lakes(i).pp(from_year)=0;
    }
    /////////////////////////////////////////////////


    /////////////////////////////////////////////////
    // Distance matrix //
    if(dist_file_name == ""){cerr << "**Use the -dist flag to provide distance matrix.\n";exit(1);}
    ifstream dist_file;
    dist_file.open(dist_file_name.c_str()); //distance_matrix.csv
    n_sources = wc_l(dist_file_name.c_str());
    if(wc_column(dist_file_name.c_str()) != n_lakes){cerr << "**Distance file should be a n_sources x n_sites table.\n";exit(1);}
    d_matrix.redim(1,n_sources,1,n_lakes);
    sources.redim(1,n_sources);

    for(int i = 1;i<=n_sources;i++)
	{
        sources(i).Dij.redim(1,n_lakes);
        sources(i).Gij.redim(1,n_lakes);
        sources(i).Xit.redim(from_year-1,to_year);
        sources(i).Xit(from_year-1)=0;
        for(int j = 1;j<=n_lakes;j++)
        {
            dist_file >> sources(i).Dij(j);
        }
    }
    dist_file.close();
    /////////////////////////////////////////////////


    /////////////////////////////////////////////////
    // Oi //
    if(Oi_file_name == ""){cerr << "**Use the -Oi flag to provide source size file.\n";exit(1);}
    if(n_sources != wc_l(Oi_file_name.c_str())) {cerr << "**Oi file must have the same number of rows as the distance matrix\n";exit(1);}
    if(wc_column(Oi_file_name.c_str()) != 1){cerr << "**Expecting only one column in Oi file.\n";exit(1);}
    ifstream o_file;
    o_file.open(Oi_file_name.c_str());

    for(int i = 1;i<=n_sources;i++)
	{
        o_file >> sources(i).Oi;
    }
    o_file.close();
    /////////////////////////////////////////////////


    /////////////////////////////////////////////////
    // init env params //
    //for(int ch=1;ch<=n_chem_var+1;ch++)
    //    chem_pars(ch)=0; //chem_pars(ch)=runif(-0.01,0.01);
    /////////////////////////////////////////////////
}

void which_sampled_or_valid()
{
   _vbc_vec<int> tmp_index_sampled(1,n_lakes);
   _vbc_vec<int> tmp_val_lakes_index(1,n_lakes);    
   _vbc_vec<int> tmp_val_inv_lakes_index(1,n_lakes); 
   int n_validation = 0;
   int n_invaded_in_validation_set = 0;
   n_sampled = 0;
   for(int i=1;i<=n_lakes;i++)
   {
      // Count validation lakes //
      if(lakes(i).status_2010 != 0)
      { 
         n_validation++;
         tmp_val_lakes_index(n_validation)=i;
         if(lakes(i).status_2010 == 2)
         {
            lakes(i).val_invaded = 1;
            n_invaded_in_validation_set++;
            tmp_val_inv_lakes_index(n_invaded_in_validation_set)=i;
         }
      }
      if( (lakes(i).discovered != 0 || lakes(i).last_abs != 0) && lakes(i).discovered != from_year) //sampled lakes (excluding the seed(s) and validation lakes)
      {
         if( (fit_pdet && lakes(i).discovered != 0) ||  !fit_pdet )
         {
             n_sampled++;
             tmp_index_sampled(n_sampled)=i;
         }
      }
   }

   cout << "Sampled "<< n_sampled << " of " << n_lakes << " lakes.\n";
   sampled_index.redim(1,n_sampled);
   for(int i = 1;i<=n_sampled;i++)
      sampled_index(i)=tmp_index_sampled(i);

   if(n_validation > 0)
   {
      cout << "Number of lakes in validation set = "<< n_validation << "\n";

      n_val_lakes = n_validation;
      val_lakes_index.redim(1,n_validation);
      for(int i = 1;i<=n_validation;i++)
         val_lakes_index(i)=tmp_val_lakes_index(i);

      val_lakes_inv_index.redim(1,n_invaded_in_validation_set);
      for(int i = 1;i<=n_invaded_in_validation_set;i++)
         val_lakes_inv_index(i)=tmp_val_inv_lakes_index(i);
   }
}

void init_state()
{
    state.redim(from_year,to_year);
    for(int t=from_year;t<=to_year;t++)
    {
        state(t).inv.redim(1,n_lakes);
        state(t).u_inv.redim(1,n_lakes);
        state(t).new_inv.redim(1,n_lakes);
        state(t).calc_pp_index.redim(1,n_lakes);
    }
}

//////////////////////////////////////////////////////////////
/////////////////// FUNCTIONS ////////////////////////////////
float grav_func(int i, int k)
{
    float tmp;
    tmp = pow(lakes(k).Wj,e_par) * pow(sources(i).Dij(k),-d_par);
    return tmp;
}
void calc_traf()
{
    //#pragma omp parallel for
    for(int i=1;i<=n_sources;i++)
    {
        sources(i).Ai=0;
        for(int k=1;k<=n_lakes;k++)
        {   
            sources(i).Gij(k) = grav_func(i,k);
            sources(i).Ai += theta_vec(k) * sources(i).Gij(k);
        }
        for(int j=1;j<=n_lakes;j++)
        {
            sources(i).Gij(j) = phi_vec(j) * sources(i).Gij(j) / (sources(i).Ai);
            if(std::isnan(sources(i).Gij(j)))
                sources(i).Gij(j) = 0;
        }
    }
}
void calc_traf_mat()
{
   //This is the bottleneck of the whole program.
   //(n_lakes*n_lakes - n_lakes)/2 * n_sources 
   // ~ 460million calcs every call.
   
   cerr << "Calculating traf mat...\n";
   //#pragma omp parallel for
    for(int i=1;i<=n_lakes;i++)
    {
        for(int j=1;j<=i-1;j++) //start at 0 since the diagonal of traf_mat is not needed.
        {
            traf_mat(i,j)=0;            
            for(int s=1;s<=n_sources;s++)
                traf_mat(i,j)+=sources(s).Gij(i)*sources(s).Gij(j)*sources(s).Oi;

            traf_mat(j,i)=traf_mat(i,j);
        }
    }
   cerr << "Finished calculating traf mat...\n";
}

// Pre-calculate an array of trafic matricies from a grid of e_par, d_par //
void calc_traf_mat_array()
{
    float step_d = (max_d - min_d)/(float(tmat_res)-1);
    float step_e = (max_e - min_e)/(float(tmat_res)-1);

    // Deal with management matrix //
    int man_dims(1);
    if(impose_management)
        man_dims = 2;
    traf_mat_array.redim(1,man_dims,1,tmat_res*tmat_res,1,n_lakes,1,n_lakes); // <- holy RAM, batman!

    cout << "Impose Management: " << impose_management << "\n";
    d_par = min_d;
    e_par = min_e;
    int m(1);
    int n(1);
    for(int di=1;di<=tmat_res;di++)
    {
        e_par = min_e;
        for(int ei=1;ei<=tmat_res;ei++)
        {
            cerr << "Pre-computing trafic matrix " << n << " of " << tmat_res*tmat_res << "\n";
            calc_traf();
            calc_traf_mat();

            for(int i=1;i<=n_lakes;i++)
            {
                for(int j=1;j<=i-1;j++) //start at 0 since the diagonal of traf_mat is not needed.
                {
                    traf_mat_array(m,n,i,j) = traf_mat(i,j);
                    traf_mat_array(m,n,j,i) = traf_mat(i,j);
                }
            }
            //cout << "\n" << d_par << " " << e_par << " " << step_d << " " << step_e << "\n";
            n++;
            e_par += step_e;
        }
        d_par += step_d;
    }
}

void calc_traf_mat_array_man()
{
    float step_d = (max_d - min_d)/(float(tmat_res)-1);
    float step_e = (max_e - min_e)/(float(tmat_res)-1);

    // Deal with management matrix //
    d_par = min_d;
    e_par = min_e;
    // -------------------------- //

    int m(2);
    int n(1);
    for(int di=1;di<=tmat_res;di++)
    {
        e_par = min_e;
        for(int ei=1;ei<=tmat_res;ei++)
        {
            cerr << "Pre-computing trafic matrix " << n << " of " << tmat_res*tmat_res << "\n";
            calc_traf();
            calc_traf_mat();

            for(int i=1;i<=n_lakes;i++)
            {
                for(int j=1;j<=i-1;j++) //start at 0 since the diagonal of traf_mat is not needed.
                {
                    traf_mat_array(m,n,i,j) = traf_mat(i,j);
                    traf_mat_array(m,n,j,i) = traf_mat(i,j);
                }
            }
            //cout << d_par << " " << e_par << " " << step_d << " " << step_e << "\n";
            n++;
            e_par += step_e;
        }
        d_par += step_d;
    } 
}

void write_traf_mat_array()
{
    ofstream tmatfile("traf_mat_array.tab");
    int n(1);
    for(int di=1;di<=tmat_res;di++)
    {
        for(int ei=1;ei<=tmat_res;ei++)
        {
            cerr << "Writing Pre-computed trafic matrix " << n << " of " << tmat_res*tmat_res << "\n";
            for(int i=1;i<=n_lakes;i++)
            {
                for(int j=1;j<=n_lakes-1;j++)
                    tmatfile << traf_mat_array(1,n,i,j) << "\t";
                tmatfile << traf_mat_array(1,n,i,n_lakes) << "\n";    
            }
            n++;
        }
    }
    tmatfile.close();
}
void read_traf_mat_array()
{
    traf_mat_array.redim(1,2,1,tmat_res*tmat_res,1,n_lakes,1,n_lakes); // <- holy RAM, batman!
    ifstream tmatfile("traf_mat_array.tab");
    int n(1);
    for(int di=1;di<=tmat_res;di++)
    {
        for(int ei=1;ei<=tmat_res;ei++)
        {
            cerr << "Reading Pre-computed trafic matrix " << n << " of " << tmat_res*tmat_res << "\n";
            for(int i=1;i<=n_lakes;i++)
            {
                for(int j=1;j<=n_lakes;j++)
                     tmatfile >> traf_mat_array(1,n,i,j);
            }
            n++;
        }
    }
    tmatfile.close();
}

// given a d_par and e_par, return the index of the nearest pre-
// calculated traffic matrix.
int lookup_traf_mat_array_index()
{
    // runs in O(1) rather than finding closest value by search.
    float step_d = (max_d - min_d)/(tmat_res-1);
    float step_e = (max_e - min_e)/(tmat_res-1);
    int ind = floor((d_par-min_d)/step_d * tmat_res) + 1;
    ind = floor((e_par-min_e)/step_e) + ind;
    return(ind);
}

void calc_traf_mat_part(_vbc_vec<int> tracked)
{
   //Calculate rows of the traf_mat only as needed.
    for(int i=1;i<=n_lakes;i++)
    {
        if(tracked(i) == 1)
        {
            for(int j=1;j<=i-1;j++)
            {
                if(tracked(j) == 1)
                {
                    //only calc if new pars or not already calc'd.
                    if(new_tr_par || traf_mat(i,j)==0)
                    {
                        traf_mat(i,j)=0;            
                        for(int s=1;s<=n_sources;s++)
                            traf_mat(i,j)+=sources(s).Gij(i)*sources(s).Gij(j)*sources(s).Oi;
                        
                        traf_mat(j,i)=traf_mat(i,j);
                    }
                    //cout << i << " " << j << ": " << traf_mat(j,i) << "\n";
                }
            }
        }
    }
}
void clear_traf_mat()
{
    for(int i=1;i<=n_lakes;i++)
    {
        for(int j=1;j<=n_lakes;j++)
        {
            traf_mat(i,j) = 0;
            traf_mat(j,i) = 0;
        }
    }
}

// Uj is the maximum propagule pressure to lake j (ie when everything is invaded).
void calc_Uj()
{
    for(int j=1;j<=n_lakes;j++)
    {
        lakes(j).Uj = 0;
        for(int i=1;i<=n_sources;i++)
        {   
            lakes(j).Uj += sources(i).Gij(j) * sources(i).Oi;
        }
    }
}
// *************************************************

void calc_pp()
{
    int lake_to_index, lake_from_index;
    for(int t=from_year+1;t<=to_year;t++)
    {
        for(int j=1;j<=state(t).n_calc_pp;j++) //n_calc_pp=the number of uninvaded+newly invaded sites
        {
            lake_to_index=state(t).calc_pp_index(j);
            lakes(lake_to_index).pp(t)=lakes(lake_to_index).pp(t-1);
            for(int i=1;i<=state(t-1).n_new_inv;i++)
            {
                lake_from_index=state(t-1).new_inv(i);
                if(traf_mat_precompute)
                    lakes(lake_to_index).pp(t) += traf_mat_array(man_dim,glb_traf_array_ind,lake_from_index,lake_to_index);
                else
                    lakes(lake_to_index).pp(t) += traf_mat(lake_from_index,lake_to_index);
            }
        }
    }
}
void calc_pp_validation(_vbc_vec<int> indicies) //for calculated pp in each relevant year for the validation lakes.
{
   int lake_to_index, lake_from_index;
   for(int l=1;l<=indicies.UBound();l++)
   {
      lake_to_index = indicies(l);
      for(int t=from_year+1;t<=to_year;t++)
      {      
         lakes(lake_to_index).pp(t) = lakes(lake_to_index).pp(t-1);
         for(int i=1;i<=state(t-1).n_new_inv;i++)
         {
            lake_from_index=state(t-1).new_inv(i);
            if(traf_mat_precompute)
                lakes(lake_to_index).pp(t) += traf_mat_array(man_dim,glb_traf_array_ind,lake_from_index,lake_to_index);
            else
                lakes(lake_to_index).pp(t) += traf_mat(lake_from_index,lake_to_index);
         }
         //cout << lake_to_index << "\t" << lakes(lake_to_index).pp(t) << "\n";
      }
   }
}

void calc_state()
{
    for(int t=from_year;t<=to_year;t++)
    {
        state(t).n_inv=0;
        state(t).n_new_inv=0;
        state(t).n_u_inv=0;
        state(t).n_calc_pp=0;
        for(int i=1;i<=n_lakes;i++)
        {
            if(t_vec(i)<=t)
            {
                if(t_vec(i)==t)
                {
                    state(t).n_new_inv++;
                    state(t).new_inv(state(t).n_new_inv)=i;
                    state(t).n_calc_pp++;
                    state(t).calc_pp_index(state(t).n_calc_pp)=i;
                }
                state(t).n_inv++;
                state(t).inv(state(t).n_inv)=i;
            }
            else
            {
                state(t).n_u_inv++;
                state(t).u_inv(state(t).n_u_inv)=i;
                state(t).n_calc_pp++;
                state(t).calc_pp_index(state(t).n_calc_pp)=i;
            }
        }
            //cout << "year\t" << t<< "\t#uninv\t" << state(t).n_u_inv << "\n";
    }   
}
void init_t()
{
    for(int i=1;i<=n_lakes;i++)
    {
        if(lakes(i).last_abs != 0 && lakes(i).invaded == 0) //observed uninvaded
        {
            t_vec(i)=runif(lakes(i).last_abs+1,to_year+2);
        }
        else if(lakes(i).last_abs != 0 && lakes(i).invaded == 1) //observed uninvaded at last_abs and observed invaded at discovered
        {
            t_vec(i)=runif(lakes(i).last_abs+1,lakes(i).discovered+1);
        }
        else if(lakes(i).discovered==from_year)
        {
            t_vec(i)=from_year;
cerr << i << "\n";
        }
        else if(lakes(i).invaded == 1) // never observed uninvaded and observed invaded at discovered
        {
             t_vec(i)=runif(from_year+1,lakes(i).discovered+1);
        }
        else
        {
            //t_vec(i)=runif(from_year+1,to_year+(to_year-from_year)); //allow for statespace in the future.
            //t_vec(i)=runif(from_year+1,to_year+2);
            t_vec(i)=2011;
        }
    }
}
void clear_t()
{
   //Set all but seed lakes to uninvaded.
   for(int i=1;i<=n_lakes;i++)
   {
      if(t_vec(i) != from_year)
         t_vec(i) = to_year + 1;
   }
}
void sample_t()
{
    int t_var=2; //+- for uniform proposal of the state variable
    for(int i=1;i<=n_lakes;i++)
    {
        if(lakes(i).last_abs != 0 && lakes(i).invaded == 0) //observed uninvaded
        {
            t_vec(i)=runif(which_max(lakes(i).last_abs+1,t_vec(i)-t_var),which_min(to_year+2,t_vec(i)+t_var+1));
        }
        else if(lakes(i).last_abs != 0 && lakes(i).invaded == 1) //observed uninvaded at last_abs and observed invaded at discovered
        {
            t_vec(i)=runif(which_max(lakes(i).last_abs+1,t_vec(i)-t_var),which_min(lakes(i).discovered+1,t_vec(i)+t_var+1) );
        }
        else if(lakes(i).invaded == 1) // never observed uninvaded and observed invaded at discovered
        {
            t_vec(i)=runif(which_max(from_year,t_vec(i)-t_var),which_min(lakes(i).discovered+1,t_vec(i)+t_var+1));
        }
        else
        {
            t_vec(i)=runif(which_max(from_year+1,t_vec(i)-t_var),which_min(to_year+2,t_vec(i)+t_var+1) ); //allow for statespace in the future.      
        }
    }
}
void sample_t_l(int i)
{
        int t_var=4; //+- for uniform proposal of the state variable
        if(lakes(i).last_abs != 0 && lakes(i).invaded == 0) //observed uninvaded
        {
            t_vec(i)=runif(which_max(lakes(i).last_abs+1,t_vec(i)-t_var),which_min(to_year+2,t_vec(i)+t_var+1));
        }
        else if(lakes(i).last_abs != 0 && lakes(i).invaded == 1) //observed uninvaded at last_abs and observed invaded at discovered
        {
            t_vec(i)=runif(which_max(lakes(i).last_abs+1,t_vec(i)-t_var),which_min(lakes(i).discovered+1,t_vec(i)+t_var+1) );
        }
        else if(lakes(i).invaded == 1) // never observed uninvaded and observed invaded at discovered
        {
            t_vec(i)=runif(which_max(from_year,t_vec(i)-t_var),which_min(lakes(i).discovered+1,t_vec(i)+t_var+1) );
        }
        else
        {
            t_vec(i)=runif(which_max(from_year+1,t_vec(i)-t_var),which_min(to_year+2,t_vec(i)+t_var+1) ); //allow for statespace in the future.      
        }
}
void sim_spread2()
{
    _vbc_vec<float> alpha(1,n_lakes);
    _vbc_vec<float> pr_I(1,n_lakes);
    int lake_index;
    bool ff;
    bool kn;

    _vbc_vec<float> unifs(1,n_lakes,from_year+1,to_year);
    _vbc_vec<int> test_lakes(1,n_lakes);
    for(int i=1;i<=n_lakes;i++)
    {
        float min_unifs = 1;
        test_lakes(i) = 0;
        alpha(i) = glb_alpha;
        pr_I(i) = calc_pr_I(i);
        for(int t=from_year+1;t<=to_year;t++)
        {
            unifs(i,t) = runif(0,1);
            if( unifs(i,t) < min_unifs ) min_unifs = unifs(i,t);
        }
        float unif = runif(0,1);
        if(pr_I(i) > unif) // || lakes(i).discovered != 0) //Invasible?
        {
            test_lakes(i) = 1;
//cout << i << " : " << pr_I(i) << " : " << unif << " : " << lakes(i).discovered << " : " << lakes(i).chem(1) << "\n";
        }
    }

    update_sim_pp(from_year,test_lakes);
    for(int t=from_year+1;t<=to_year;t++)
    {
        state(t).n_inv=state(t-1).n_inv;
        state(t).n_u_inv=0;
        state(t).n_new_inv=0;
        for(int i=1;i<=state(t-1).n_u_inv;i++)
        {
            lake_index=state(t-1).u_inv(i);
            //update_pp_l_hood(lake_index,t); // TEST to determine if equiv to update_sim_pp.
            if(test_lakes(lake_index) == 1)
            {
                if(force_invasions && lakes(lake_index).discovered==t && lakes(lake_index).last_abs < t)
                    ff = TRUE;
                else
                    ff = FALSE;
                if( 1-exp(- pow(alpha(lake_index) * lakes(lake_index).pp(t) +
                         gamma_par, c_par ) )  > unifs(lake_index,t) || ff)
                {
                    state(t).n_inv++;
                    state(t).n_new_inv++;
                    state(t).new_inv(state(t).n_new_inv)=lake_index;
                    state(t).inv(state(t).n_inv)=lake_index;
                    t_vec(lake_index)=t;
                    
                } else {
                    state(t).n_u_inv++;
                    state(t).u_inv(state(t).n_u_inv)=lake_index;
                    t_vec(lake_index)=to_year+1;
                }
            }else{
                state(t).n_u_inv++;
                state(t).u_inv(state(t).n_u_inv)=lake_index;
                t_vec(lake_index)=to_year+1;
            }
        }     
        if(t<to_year)
            update_sim_pp(t,test_lakes);
    }
}

void sim_spread()
{
    _vbc_vec<float> alpha(1,n_lakes);
    int lake_index;
    bool ff;
    bool kn;
    // determine which lakes we need to fill in traf_mat for:
    _vbc_vec<float> unifs(1,n_lakes,from_year+1,to_year);
    _vbc_vec<int> test_lakes(1,n_lakes);
    for(int i=1;i<=n_lakes;i++)
    {
        float min_unifs = 1;
        test_lakes(i) = 0;
        alpha(i)=calc_alpha(i);
        for(int t=from_year+1;t<=to_year;t++)
        {
            unifs(i,t) = runif(0,1);
            if( unifs(i,t) < min_unifs ) min_unifs = unifs(i,t);
        }
        if( 1-exp(- pow(alpha(i) * lakes(i).Uj + gamma_par, c_par ) ) >  \
            min_unifs || lakes(i).discovered != 0 || lakes(i).last_abs != 0 )
        {
            //cout << lakes(i).last_abs << "++++\n";
            test_lakes(i) = 1;
        }
        test_lakes(i) = 1; // TTMP
    }
    //if(!one_time_traf)
    //   calc_traf_mat_part(test_lakes);
            //cout << "------------------------------------\n";
    man_dim=1;
    update_sim_pp(from_year,test_lakes);
    for(int t=from_year+1;t<=to_year;t++)
    {
        state(t).n_inv=state(t-1).n_inv;
        state(t).n_u_inv=0;
        state(t).n_new_inv=0;
        for(int i=1;i<=state(t-1).n_u_inv;i++)
        {
            lake_index=state(t-1).u_inv(i);
            //update_pp_l_hood(lake_index,t); // TEST to determine if equiv to update_sim_pp.
            if(test_lakes(lake_index) == 1)
            {
                //cout << "*********" << lakes(lake_index).pp(t) << "\n";
                // invade stochastically and restrict to observed pattern
                if(force_invasions && lakes(lake_index).discovered==t && lakes(lake_index).last_abs < t)
                    ff = TRUE;
                else
                    ff = FALSE;
                if(lakes(lake_index).discovered == 0 && only_known)
                    kn = FALSE;
                else
                    kn = TRUE;

                if(  (( 1-exp(- pow(alpha(lake_index) * lakes(lake_index).pp(t) +
                         gamma_par, c_par ) )  > unifs(lake_index,t) ) \
                         || ff ) &&
                         kn) //only let known sites invade
                {
                    state(t).n_inv++;
                    state(t).n_new_inv++;
                    state(t).new_inv(state(t).n_new_inv)=lake_index;
                    state(t).inv(state(t).n_inv)=lake_index;
                    t_vec(lake_index)=t;
                    
                } else {
                    state(t).n_u_inv++;
                    state(t).u_inv(state(t).n_u_inv)=lake_index;
                    t_vec(lake_index)=to_year+1;
                }
            }else{
                state(t).n_u_inv++;
                state(t).u_inv(state(t).n_u_inv)=lake_index;
                t_vec(lake_index)=to_year+1;
            }
        }
        ////////////////////////////////////////////////////////
        // Undoing invasions to compensate for forcing lakes: //
        ////////////////////////////////////////////////////////
        if(compensate_force)
        {
            int n_forced(0),n_not_forced(0);
            _vbc_vec<int> forced(1,state(t).n_new_inv);
            _vbc_vec<int> not_forced(1,state(t).n_new_inv);
            _vbc_vec<int> not_forced_lake_index(1,state(t).n_new_inv);
            _vbc_vec<float> PInv(1,state(t).n_new_inv);
            for(int i=1;i<=state(t).n_new_inv;i++)
            {
                lake_index = state(t).new_inv(i);
                PInv(i) = 1-exp(- pow(alpha(lake_index) *lakes(lake_index).pp(t) + \
                        gamma_par, c_par ));
                if(lakes(lake_index).discovered == t && \
                    PInv(i) < unifs(lake_index,t))
                {
                    n_forced++;
                    forced(n_forced)=i;
                } else {
                    n_not_forced++;
                    not_forced_lake_index(n_not_forced) = lake_index;
                    not_forced(n_not_forced)=i;
                }
            }
            cout << "\n----------------------------------------------------\n" \
               << t << " " << n_forced << " forced invasions. " << n_not_forced << "unforced invasions. Difference: " << n_not_forced-n_forced << "\n";
            
            for(int f=1;f<=n_forced;f++)
            {
                if(n_not_forced > 0)
                {
                    float min_dist(1000), dist, min_rel_dist;
                    int which_min;
                    for(int k=1;k<=n_not_forced;k++)
                    {
                        //find closest PInv
                        dist = pow(PInv(forced(f))-PInv(not_forced(k)),2);
                        if(dist < min_dist)
                        {
                            min_dist = dist;
                            min_rel_dist = PInv(forced(f))-PInv(not_forced(k));
                            which_min = k; 
                        }
                    }
                    // Swap //
                    //cout <<  f  << " swapped for " << not_forced(which_min);
                    //cout << ". Diff PInv: " << PInv(forced(f)) << " - " << PInv(not_forced(which_min));
                    //cout << " = " << min_rel_dist << "\n"; 

                    //for(int jj=1;jj<=n_not_forced;jj++)
                    //    cout << not_forced_lake_index(jj) << " :: ";
                    //cout << "\n";
                    // add to uninv //
                    state(t).n_u_inv++;
                    state(t).u_inv(state(t).n_u_inv)=not_forced_lake_index(which_min);
                    t_vec(not_forced_lake_index(which_min))=to_year+1;

                    // remove from inv //
                    state(t).new_inv( not_forced(which_min) ) = \
                        state(t).new_inv(state(t).n_new_inv);
                    state(t).inv( not_forced(which_min) ) = state(t).inv(state(t).n_inv);
                    state(t).n_inv--;
                    state(t).n_new_inv--;

                    // remove from not_forced list.
                    PInv(which_min) = PInv(n_not_forced);
                    not_forced(which_min) = not_forced(n_not_forced);
                    not_forced_lake_index(which_min) = not_forced_lake_index(n_not_forced);
                    n_not_forced--;                    
                }
            }
            
        }
        ////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////
        //toggle management on and off
        if(impose_management)
        {
            if(t < 2013)
                man_dim=1;
            else
                man_dim=2;

            // update the full pp timeline so as not to have residual extra pp //
            if(t == 2013)
            {
                calc_state();
                calc_pp();
            }
        }
        if(t<to_year)
            update_sim_pp(t,test_lakes);
      //cout << state(t-1).n_inv << " ";
    }
    //cout << fit_pdet << "::" << compensate_force << "::" << one_time_traf << "\n";
    //cout << "\r" << to_year << "=====" << float(state(to_year).n_inv)/float(n_lakes);
}
void update_sim_pp(int t,_vbc_vec<int> test_lakes)
{
    //add pp from each newly invaded lake at time t to each uninvaded lake
    int to_lake_index;
    int from_lake_index;
    //could just do this over those lakes that we will either test or need in likelihood.
    for(int i=1;i<=state(t).n_u_inv;i++)
    {
        to_lake_index=state(t).u_inv(i);
        if(test_lakes(to_lake_index) == 1)
        {   
           lakes(to_lake_index).pp(t+1) = lakes(to_lake_index).pp(t);
           for(int j=1;j<=state(t).n_new_inv;j++)
           {
               from_lake_index=state(t).new_inv(j);
               if(traf_mat_precompute)
                    lakes(to_lake_index).pp(t+1) += traf_mat_array(man_dim,glb_traf_array_ind,from_lake_index,to_lake_index);
               else
                    lakes(to_lake_index).pp(t+1) += traf_mat(from_lake_index,to_lake_index);
           }
        }
    }
}
void update_pp_l_hood(int lake_index,int t) //called from l_hood to calc pp where it is missing from sim
{
    int to_lake_index = lake_index;
    int from_lake_index;
    lakes(to_lake_index).pp(t) = lakes(to_lake_index).pp(t-1);
    //cout << to_lake_index << "\t" << t << "\t" << lakes(to_lake_index).pp(t) << "\n";
    for(int j=1;j<=state(t-1).n_new_inv;j++)
    {
        from_lake_index = state(t-1).new_inv(j);
        if(from_lake_index != to_lake_index){ //should not receive pp from itself!
            if(traf_mat_precompute)
                lakes(to_lake_index).pp(t) += traf_mat_array(man_dim,glb_traf_array_ind,from_lake_index,to_lake_index);
            else
                lakes(to_lake_index).pp(t) += traf_mat(from_lake_index,to_lake_index);
            //cout << traf_mat(from_lake_index,to_lake_index) << " @ ";
        }
    }
   //cout << "\n";
}


float l_hood()
{
    float lh(0); float tmp_lh;
    int lake_index;
    for(int i=1;i<=n_sampled;i++)
    {
        lake_index=sampled_index(i);
        float alpha=calc_alpha(lake_index); //insert this inside lake loop.    
        if(lakes(lake_index).last_abs != 0) // observed uninvaded
        {
            for(int t=from_year+1;t<=lakes(lake_index).last_abs;t++)
            {
                if(t_vec(lake_index) < t)
                  update_pp_l_hood(lake_index,t);
// lake is simmed invaded later than where the pp record goes to //
// How could this happen? PP should be calced for all lakes that.//

                //Prob uninv to that year
                lh -= pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par); 
            }
            if(lakes(lake_index).discovered != 0) // discovered invaded
            {
                tmp_lh=0;
                for(int t=lakes(lake_index).last_abs+1;t<=lakes(lake_index).discovered;t++)
                {
                   if(t_vec(lake_index)<t)
                       update_pp_l_hood(lake_index,t);

                    tmp_lh -= pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par);
                }
                lh+= log(1-exp(tmp_lh));
            }
        }else{
            tmp_lh=0;
            //never observed uninvaded, but discovered at T
            for(int t=from_year+1;t<=lakes(lake_index).discovered;t++)
            {
               if(t_vec(lake_index)<t)
                  update_pp_l_hood(lake_index,t);

               tmp_lh -= pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par);
            }
            lh += log(1-exp(tmp_lh));
        }
 
/*    if(lake_index ==134) //pp only getting calc'd in 1990? wtf?
    {
        cout << lake_index << ": " << lakes(lake_index).discovered << ":: " << lakes(lake_index).last_abs << ":: "<< lh << "\n";
        cout << lakes(lake_index).pp(lakes(lake_index).discovered) << " :: " << lakes(lake_index).pp(lakes(lake_index).last_abs-2) << "\n";
        cout << t_vec(134) << "***\n";
        for(int i=1;i<=state(2009).n_u_inv;i++)
            cout << state(2009).u_inv(i) << "---\n";
        for(int i=1;i<=state(2009).n_inv;i++)
            cout << state(2009).inv(i) << "+++\n";

    }
*/
    }
    return lh;
}

// Likelihood for the pressence only model //
float l_hood_detp()
{
    float lh(0); float tmp_lh;
    int lake_index;
    float alpha;
    float pdet_;
    //for(int i=1;i<=n_sampled;i++)
    //{
      //lake_index=sampled_index(i);
    for(int i=1;i<=n_lakes;i++)
    {
         lake_index=i;
         alpha = calc_alpha(lake_index);
         pdet_ = calc_pdet(lake_index);
         tmp_lh=0;
         if((!include_undetected_sites && lakes(lake_index).discovered == 0) || lakes(lake_index).discovered == from_year)
               continue;

         //never observed uninvaded, but discovered at T
         for(int t=from_year+1;t<=to_year;t++) 
         {
            //cout << i << " " << t << " " << t_vec(i) << " " << lakes(i).discovered << " " << tmp_lh << "\n";
            //if(t_vec(lake_index) < t)   
            //   update_pp_l_hood(lake_index,t);
            // The strategy here is that we use (1-poe) 
            // for t>t_vec(lake_index), poe for t = t_vec(lake_index) 
            // and (1-pdet) from t=t_vec(lake_index) to t<lakes(lake_index).discovered, 
            // pdet for t=lakes(lake_index).discovered.
            // 
            update_pp_l_hood(lake_index,t);

            //Prob uninv that year
            if(t < t_vec(lake_index))
            {               
               tmp_lh += - pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par); 
            }

            //Prob inv IN that year
            if(t == t_vec(lake_index))
            {  
               tmp_lh += log(1 - exp( -pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par) ));
            }

            //Prob of not detecting
            if(t >= t_vec(lake_index) && //past the time of invasion
               t != lakes(lake_index).discovered && //not the discovery year
               lakes(lake_index).discovered != from_year) //not a seed lake
            {
                 tmp_lh += log(1 - pdet_); 
            }

            //Prob of detecting
            if(t == lakes(lake_index).discovered) 
            {
                //To handle the non-forced simulations//
                if(t_vec(lake_index) > lakes(lake_index).discovered)
                {
                    //remove the Prob uninv that year
                    tmp_lh += pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par); 
                    // add Prob inv
                    tmp_lh += log(1 - exp( -pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par) ));
                }
                tmp_lh += log(pdet_); //Prob of detecting in the year of detection.
                break;
            }

         }
         lh += tmp_lh;   
    }
//cout << "\n---------------------------------------------------------------\n";
    return lh;
}
float l_hood_detp2() //integrate over time (using mean field prob inv)
{
    float lh(0); float tmp_lh; float tmp_lh_running;
    int lake_index;
    int len;
    float alpha;
    //pdet = gen_pdet;
    for(int i=1;i<=n_lakes;i++)
    {
        lake_index = i;
        tmp_lh = 1;
        int to = to_year;
        tmp_lh_running = 0;
        if(lakes(lake_index).discovered != from_year)
        {
            if(lakes(lake_index).discovered != 0)
            {   
                for(int t_inv=from_year+1;t_inv<=lakes(lake_index).discovered;t_inv++)
                {
                    //int t_inv = t_vec(lake_index);
                    tmp_lh = 1;
                    for(int t=from_year+1;t<=lakes(lake_index).discovered;t++)
                    {
                        if(t < t_inv)
                            tmp_lh *= 1 - lakes(lake_index).PInv_inst(t);
                        if(t == t_inv || t == lakes(lake_index).discovered)
                            tmp_lh *= lakes(lake_index).PInv_inst(t); 
                        if(t >= t_inv && t != lakes(lake_index).discovered)
                            tmp_lh *= 1 - pdet;
                    }
                    tmp_lh *= pdet;
                    tmp_lh_running += tmp_lh;
                }
            } else {
                for(int t_inv=from_year+1;t_inv<=to_year;t_inv++)
                {
                    //int t_inv = t_vec(lake_index);
                    tmp_lh = 1;
                    for(int t=from_year+1;t<=to_year;t++)
                    {
                        if(t < t_inv)
                            tmp_lh *= 1 - lakes(lake_index).PInv_inst(t);
                        if(t == t_inv)
                            tmp_lh *= lakes(lake_index).PInv_inst(t); 
                        if(t >= t_inv)
                            tmp_lh *= 1 - pdet;
                    }
                    tmp_lh_running += tmp_lh;
                 }
            }
            lh += log(tmp_lh_running);
        }
    }
    return lh;
}
float l_hood_detp3() //Integrate over time (with simulated states).
{
    float lh(0); float tmp_lh; float tmp_lh_running; float tmp;
    int lake_index;
    int len;
    float alpha;
    float pdet_;
    for(int i=1;i<=n_lakes;i++)
    {
        lake_index = i;
        tmp_lh_running = 0;
        alpha = calc_alpha(lake_index);
        pdet_ = calc_pdet(lake_index);
        if(lakes(lake_index).discovered != from_year)
        {
            // Observed (detections)
            if(lakes(lake_index).discovered != 0)
            {   
                for(int t_inv=from_year+1;t_inv<=lakes(lake_index).discovered;t_inv++)
                {
                    //int t_inv = t_vec(lake_index);
                    tmp_lh = 1; // The summed probability across states //
                    for(int t=from_year+1;t<=lakes(lake_index).discovered;t++)
                    {
                        //update propagule pressure //
                        if(t_inv == from_year+1)
                            update_pp_l_hood(lake_index,t);
                        //////////////////////////////

                        if(t < t_inv)
                            tmp_lh *= exp( -pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par) ); //log(1-p_inv)
                        if(t == t_inv)
                            tmp_lh *= 1 - exp( -pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par) );//lakes(lake_index).PInv_inst(t); 
                        if(t >= t_inv && t != lakes(lake_index).discovered)
                            tmp_lh *= 1 - pdet_;
                    }
                    tmp_lh *= pdet_;

                    //cout << lake_index << " " << t_inv << " " << pdet << " " <<1 - exp( -pow(alpha*lakes(lake_index).pp(to_year)+gamma_par,c_par) ) << " " << tmp_lh << "\n";
                    tmp_lh_running += tmp_lh;
                }
                    tmp_lh_running = tmp_lh_running; // (lakes(lake_index).discovered - from_year + 1);
                //cout << "------------\n" << tmp_lh_running << "\n---------------\n";
            } else {
            // All sites with no detections 
                for(int t_inv=from_year+1; t_inv<=to_year+1; t_inv++)
                {
                    //int t_inv = t_vec(lake_index);
                    tmp_lh = 1;
                    for(int t=from_year+1;t<=to_year;t++)
                    {
                        //update propagule pressure //
                        if(t_inv == from_year+1)
                            update_pp_l_hood(lake_index,t);
                        //////////////////////////////

                        if(t < t_inv)
                            tmp_lh *= exp( -pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par) ); //log(1-p_inv)
                        if(t == t_inv)
                            tmp_lh *= 1 - exp( -pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par) );
                        if(t >= t_inv)
                            tmp_lh *= 1 - pdet_;
                    }
                    //cout << lake_index << " " << t_inv << " " << pdet << " " <<1 - exp( -pow(alpha*lakes(lake_index).pp(to_year)+gamma_par,c_par) ) << " " << tmp_lh << "\n";
                    tmp_lh_running += tmp_lh;
                }
                tmp_lh_running = tmp_lh_running; // (to_year + 1 - from_year + 1);
                //cout << "------------\n" << tmp_lh_running << "\n---------------\n";
            }
            //cout << i << " " << tmp_lh_running << " " << lh << "\n";
            lh += log(tmp_lh_running);
        }
    }
    return lh;
}

float l_hood_detp4() //Integrate over time (with simulated states) after Leung & Mandrak.
{
    float lh(0); float tmp_lh; float tmp_lh_running; float tmp;
    int lake_index;
    int len;
    float alpha = glb_alpha;
    float pr_I;
    float pdet_;
    for(int i=1;i<=n_lakes;i++)
    {
        lake_index = i;
        tmp_lh_running = 0;
        pr_I = calc_pr_I(lake_index);
        pdet_ = calc_pdet(lake_index);
        if(lakes(lake_index).discovered != from_year)
        {
            // Observed (detections)
            if(lakes(lake_index).discovered != 0)
            {   
                for(int t_inv=from_year+1;t_inv<=lakes(lake_index).discovered;t_inv++)
                {
                    tmp_lh = 1; // The summed probability across states //
                    for(int t=from_year+1;t<=lakes(lake_index).discovered;t++)
                    {
                        //update propagule pressure //
                        if(t_inv == from_year+1)
                            update_pp_l_hood(lake_index,t);
                        //////////////////////////////

                        if(t < t_inv)
                            tmp_lh *= exp( -pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par) ); //log(1-p_inv)
                        if(t == t_inv)
                            tmp_lh *= 1 - exp( -pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par) );//lakes(lake_index).PInv_inst(t); 
                        if(t >= t_inv && t != lakes(lake_index).discovered)
                            tmp_lh *= 1 - pdet_;
                    }
                    tmp_lh *= pdet_;
                    tmp_lh *= pr_I;
                    tmp_lh_running += tmp_lh;
                }
            } else {
            // All sites with no detections 
                for(int t_inv=from_year+1; t_inv<=to_year+1; t_inv++)
                {
                    //int t_inv = t_vec(lake_index);
                    tmp_lh = 1;
                    for(int t=from_year+1;t<=to_year;t++)
                    {
                        //update propagule pressure //
                        if(t_inv == from_year+1)
                            update_pp_l_hood(lake_index,t);
                        //////////////////////////////

                        if(t < t_inv)
                            tmp_lh *= exp( -pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par) ); //log(1-p_inv)
                        if(t == t_inv)
                            tmp_lh *= 1 - exp( -pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par) );
                        if(t >= t_inv)
                            tmp_lh *= 1 - pdet_;
                    }
                    if(t_inv > to_year)
                    {
                        tmp_lh *= (1-pr_I);
                    } else {
                        tmp_lh *= pr_I;
                    }
                    tmp_lh_running += tmp_lh;
                    //cout << lake_index << " : "<< t_inv << " : "<< tmp_lh << " : " << pr_I << " : "<<  tmp_lh_running <<"\n";
                }
            }
            lh += log(tmp_lh_running);
        }
    }
    return lh;
}

float l_hood_detp5() //after Leung & Mandrak.
{
    float lh(0); float tmp_lh; float tmp_lh_running; float tmp;
    int lake_index;
    int len;
    float alpha = glb_alpha;
    float pr_I;
    float pdet_;
    for(int i=1;i<=n_lakes;i++)
    {
        lake_index = i;
        tmp_lh_running = 0;
        pr_I = calc_pr_I(lake_index);
        pdet_ = calc_pdet(lake_index);
        if(lakes(lake_index).discovered != from_year)
        {
            // Observed (detections)
            if(lakes(lake_index).discovered != 0)
            {   
                    int t_inv = t_vec(lake_index);
                    tmp_lh = 1; // The summed probability across states //
                    for(int t=from_year+1;t<=lakes(lake_index).discovered;t++)
                    {
                        //update propagule pressure //
                        if(t_inv == from_year+1)
                            update_pp_l_hood(lake_index,t);
                        //////////////////////////////

                        if(t < t_inv)
                            tmp_lh += -pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par); //log(1-p_inv)
                        if(t == t_inv)
                            tmp_lh += log(1 - exp( -pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par) ));//lakes(lake_index).PInv_inst(t); 
                        if(t >= t_inv && t != lakes(lake_index).discovered)
                            tmp_lh += log(1 - pdet_);
                    }
                    tmp_lh += log(pdet_);
                    tmp_lh += log(pr_I);
                    tmp_lh_running = tmp_lh;
            } else {
            // All sites with no detections 
                    int t_inv = t_vec(lake_index);
                    //int t_inv = t_vec(lake_index);
                    tmp_lh = 0;
                    for(int t=from_year+1;t<=to_year;t++)
                    {
                        //update propagule pressure //
                        if(t_inv == from_year+1)
                            update_pp_l_hood(lake_index,t);
                        //////////////////////////////

                        if(t < t_inv)
                            tmp_lh += -pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par); //log(1-p_inv)
                        if(t == t_inv)
                            tmp_lh += log( 1 - exp( -pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par) ));
                        if(t >= t_inv)
                            tmp_lh += log(1 - pdet_);
                    }
                    if(t_inv > to_year)
                    {
                        tmp_lh += log(1-pr_I);
                    } else {
                        tmp_lh += log(pr_I);
                    }
                    tmp_lh_running = tmp_lh;
                    //cout << lake_index << " : "<< t_inv << " : "<< tmp_lh << " : " << pr_I << " : "<<  tmp_lh_running <<"\n";
            }
            lh += tmp_lh_running;
        }
    }
    return lh;
}


void mean_field()
{
    // PInv is the cummulative probability of a lake becomming invaded //
    // PInv_inst is the instantaneaous probability in a given year //
    _vbc_vec<float> alpha(1,n_lakes);
    float tmp_pp;
    for(int i=1;i<=n_lakes;i++)
    {
        //seed lakes always have PInv = 1
        if(lakes(i).discovered == from_year)
            lakes(i).PInv(from_year) = 1;
        else
            lakes(i).PInv(from_year) = 0;  
    }

    for(int t=from_year+1;t<=to_year;t++)
    {
        for(int i=1;i<=n_lakes;i++)
        {
            tmp_pp = 0;
            alpha(i)=calc_alpha(i);
            for(int j=1;j<=n_lakes;j++)
            {
                if(i != j)
                    tmp_pp += traf_mat(i,j) * lakes(j).PInv(t-1);
            }

            lakes(i).PInv_inst(t) = 1 - exp(- pow(alpha(i) * tmp_pp + gamma_par, c_par ));
            //cout << i << " " << t << " " << lakes(i).PInv_inst(t) <<"\n";
            lakes(i).PInv(t) = lakes(i).PInv(t-1) + (1-lakes(i).PInv(t-1))*lakes(i).PInv_inst(t);

            // Incorporate data //
            /*
            if(lakes(i).discovered == t)
                lakes(i).PInv(t) = 1; //Known to be invaded
            if(lakes(i).last_abs >= t)
                lakes(i).PInv(t) = 0; //Known to be uninvaded

            */
        }        
    }
    // Incorporates last_abs and discovered information //
    // for last_abs, hold PInv = 0 for all time before last_abs.
    // for discovered, set PInv = 1 for that time for the purpose of calculating everything else,
    // but retain the unconstrained PInv for the likelihood calculation (same for last_abs).
}

float calc_pr_I(int i)
{
    float z=chem_pars(1); //intercept
    if(est_env > 0)
    {
       for(int ch=2;ch<=est_env+1;ch++)
       {
           z += chem_pars(ch)*lakes(i).chem(ch-1);
       }
       if(alpha_order==2) // 2nd order polynomial //
       {
           for(int ch=2;ch<=est_env+1;ch++)
           {
               z += chem_pars(est_env+ch)*pow(lakes(i).chem(ch-1),2);
           }            
       }
    }
    float pr_I = (1/(1+exp(-z))); //Simple Logistic
    return pr_I;
}
float calc_alpha(int i)
{
   if(env)
   {
       float z=chem_pars(1); //intercept
       if(est_env > 0)
       {
           for(int ch=2;ch<=est_env+1;ch++)
           {
               z += chem_pars(ch)*lakes(i).chem(ch-1);
           }
           if(alpha_order==2) // 2nd order polynomial //
           {
               for(int ch=2;ch<=est_env+1;ch++)
               {
                   z += chem_pars(est_env+ch)*pow(lakes(i).chem(ch-1),2);
               }            
           }
       }
       float alpha = -log(1-(1/(1+exp(-z)))); //see notebook for derivation of functional form
       //float alpha = z;
       //cout << "lake "<<i << "\t" << alpha <<"\n";
       return alpha;
   }else{
      return glb_alpha;
   }
}
float calc_pdet(int i)
{
   if(hetero_pdet)
   {
       float z=pdet_pars(1); //intercept
       if(est_pdet > 0)
       {
           for(int ch=2;ch<=est_pdet+1;ch++)
           {
               z += pdet_pars(ch)*lakes(i).pdet_var(ch-1);
           }
       }
       float pdet_ = 1/(1+exp(-z)); //logistic

       return pdet_;
   }else{
      return pdet;
   }
}


// Gives the probability that lake_index is invaded by year but not yet detected.
float pr_inv_not_det(int lake_index,int year)
{
    float tmp_lh;
    float pr_inv_undet=0;
    float alpha = calc_alpha(lake_index);
    float pdet_ = calc_pdet(lake_index);
    for(int t_inv=from_year+1; t_inv<=year; t_inv++)
    {
        tmp_lh = 1;
        for(int t=from_year+1;t<=t_inv;t++)
        {
            //update propagule pressure //
            update_pp_l_hood(lake_index,t);
            //////////////////////////////

            if(t < t_inv)
                tmp_lh *= exp( -pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par) ); //log(1-p_inv)
            if(t == t_inv)
                tmp_lh *= 1 - exp( -pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par) );
        }
        tmp_lh *= pow(1 - pdet_,year-t_inv+1);
        pr_inv_undet += tmp_lh;
    }
    return(pr_inv_undet);
}


/// Wrapper for l_hood that sims spread n_sim times
/// and to smooth a quasi-likelihood function for fitting
/// MLE estimates.
float MLE_l_hood(_vbc_vec<float> * pars, _vbc_vec<float> * dat)
{
   _vbc_vec<float> tmplhood(1,n_sim_for_smooth);
   _vbc_vec<float> params = *pars;
   link_pars(parnames,(*pars));
   // Parameter bounds //
   if(d_par <=0 || e_par < 0 || c_par <= 0 || gamma_par < 0 || glb_alpha < 0 || pdet < 0 || pdet > 1)
      return(10000000);

   if(fixed_d == 0)
   {
      if(!one_time_traf)
      {
         calc_traf();
         calc_Uj();
         new_tr_par=TRUE;
         //clear_traf_mat();
      }
      //calc_traf_mat();
      //calc_pp();
   }else{
      //*pars(1)=fixed_d;
      d_par=fixed_d;
   }
   
   //cout << l_hood_detp() << "\n";
   for(int i=1;i<=n_sim_for_smooth;i++)
   {
        sim_spread();
        if(fit_pdet)
            tmplhood(i) = l_hood_detp();
        else
            tmplhood(i) = l_hood();
        //cout << i << "\t" << tmplhood(i) <<"\n";
        new_tr_par=FALSE;
   }

   float qll=average(tmplhood);
   //qll += dbeta(pdet,30,1,1);
   //cout << dbeta(pdet,30,1,1) << " --- ";

   for(int i=1;i<=params.UBound();i++)
      cout << setw(15) << params(i) <<"\t";
   cout << "===== (" << float(state(to_year).n_inv)/float(n_lakes) << ") ";
   if(sim)
      cout << "[" << gen_prop << "] ";
   cout << qll <<"\n";
   if(std::isnan(qll))
      return(10000000);


    // Print out distribution of alpha values at MLE
    /*
    ofstream alphas_file;
    alphas_file.open("output/alphas.tab",std::fstream::app);
    for(int i=1;i<=n_lakes;i++)
        alphas_file << calc_alpha(i) << "\t";
    alphas_file << "\n";
    alphas_file.close();
    */

   return(-qll);//-tve because simplex is a minimizer
}


/// MCMC LIB ////
void likelihood_wrapperMCMC(_vbc_vec<float> * params, float * l,int dim)
{
	*l = - MLE_l_hood(params, params);
}
void likelihood_wrapperMCMC_MD(_vbc_vec<float> * pars, float * l,int dim)
{
   float llmd;

   link_pars(parnames,(*pars));
   //calc_pp();
   //mean_field(); <- Good idea, but has biases for some network configurations.

   n_sim_for_smooth = 1;
   _vbc_vec<float> tmplhood(1,n_sim_for_smooth);
   for(int i=1;i<=n_sim_for_smooth;i++)
   {
        sim_spread();
        if(fit_pdet)
            tmplhood(i) = l_hood_detp3();
        else
            tmplhood(i) = l_hood();
        //cout << i << "\t" << tmplhood(i) <<"\n";
   }
   cout << "== (" << float(state(to_year).n_inv)/float(n_lakes) << ") ";
   if(sim)
      cout << "[" << gen_prop << "] ";

   llmd = average(tmplhood);// - log(glb_alpha);

   int n_par = (*pars).UBound();

   for(int i=1;i<=n_par;i++)
      cout << setw(8) << (*pars)(i) <<" ";

   cout << llmd <<"\n";

   *l = llmd;
}

//runs the join model after Leung and Mandrak 2007.
void likelihood_wrapperMCMC_MD_joint(_vbc_vec<float> * pars, float * l,int dim)
{
   float llmd;
   link_pars(parnames,(*pars));

   n_sim_for_smooth = 1;

   sim_spread2();
   llmd = l_hood_detp4();

   cout << "== (" << float(state(to_year).n_inv)/float(n_lakes) << ") ";
   if(sim)
      cout << "[" << gen_prop << "] ";

   int n_par = (*pars).UBound();

   for(int i=1;i<=n_par;i++)
      cout << setw(8) << (*pars)(i) <<" ";

   cout << llmd <<"\n";

   *l = llmd;
}
bool restrict_MCMC(float param,int dim)
{
   //link_pars(parnames,param);
   //if(d_par <=0 || e_par < 0 || c_par <= 0 || gamma_par < 0 || glb_alpha < 0 || pdet < 0 || pdet > 1)
   //   return TRUE;

   return FALSE;
}
float prior_MD(_vbc_vec<float> x, int dim)
{
   //if(!env)
   //   return -log(glb_alpha); //prior on alpha
   //else
	   return 0; //uninformative prior (log)
}
bool restrict_MCMC_MD(_vbc_vec<float> param)
{
   link_pars(parnames,param);
   if(d_par <=0 || e_par < 0 || c_par <= 0 || gamma_par < 0 || glb_alpha < 0 || pdet < 0 || pdet > 1) 
      return TRUE;

   if(traf_mat_precompute){
       if(d_par < min_d || e_par < min_e || d_par > max_d || e_par > max_e) 
          return TRUE;
   }

   return FALSE;
}
float prior(float x, int dim)
{
	return 0; //uninformative prior (log)
}



///////////////// Posterior Spead Sim ///////////////
void sim_spread_posterior()
{
   ifstream post_file("output/thinned.mcmc");
   int index;
   for(int i=1;i<=post_length;i++)
   {
      post_file >> index;
      post_file >> d_par;
      post_file >> e_par;      
      post_file >> c_par;
      if(env)
      {
         for(int i=1;i<=n_chem_var+1;i++)
            post_file >> chem_pars(i);
      }
      calc_traf();
      calc_traf_mat();

      sim_spread();
      for(int i=1;i<=10;i++)
      {
         sim_spread();
         write_t();
      }
      cout << i << " of "<< post_length<< "\n";
   }

   post_file.close();
   t_file.close();
}


/////////////// Prediction probability vectors //////////////

_vbc_vec<float> predict_p(_vbc_vec<float> params,_vbc_vec<int> indicies,int m_pars) //pars,indicies
{
   _vbc_vec<float> pred_p(1,m_pars,1,n_val_lakes);
   ofstream val_sim_file;
   val_sim_file.open("output/val_sim_props.tab");

   for(int m=1;m<=m_pars;m++)
   {
      cout << m << " of " << m_pars << "\n";
      //for(int i=1;i<=n_pars,)
      //link_pars(parnames,(*parval));
      d_par = params(m,1);
      e_par = params(m,2);
      //e_par = 1;
      c_par = params(m,3);
      //gamma_par = params(m,4);
      //glb_alpha = params(m,5);
      calc_traf();
      clear_traf_mat();
         //calc_traf_mat();
         //calc_pp();
      
      if(env)
      {
         for(int i=1;i<=n_chem_var+1;i++){
            chem_pars(i)=params(m,3+i);
            cout << chem_pars(i) << "\n";}
      }


      // -- Using the simulation method -- //
      // prob of inv is calculated as the proportion of times invaded via simulation // 
      _vbc_vec<float> prop_val_invaded(1,n_val_lakes);
      for(int i=1;i<=n_val_lakes;i++)
         prop_val_invaded(i) = 0;

      int n_sims = 1000;
      for(int s=1;s<=n_sims;s++)
      {
         
         sim_spread();
         new_tr_par = FALSE;
         for(int i=1;i<=n_val_lakes;i++)
         {
            if(t_vec(indicies(i)) != to_year+1)
               prop_val_invaded(i) += 1;
         }
      }
      cout << "\n";
      for(int i=1;i<=n_val_lakes;i++)
      {
            prop_val_invaded(i) = float(prop_val_invaded(i))/float(n_sims);
            val_sim_file << prop_val_invaded(i) << "\t";
            cout << prop_val_invaded(i) << "\t";
      }
      cout << "\n";
      val_sim_file << "\n";
      // --------------------------------------------- //


      calc_pp_validation(indicies); // fill in all pp values for validation lakes
                                    // since calc_pp() is optimized to only calc pp for uninvaded lakes.
      int cal_from;
      for(int i=1;i<=n_val_lakes;i++)
      {
         float alpha=calc_alpha(indicies(i)); //insert this inside lake loop.
         // 1 - prob of remaining uninvaded during the period since last observed absence //
         if(lakes(indicies(i)).last_abs != 0)
            cal_from = lakes(indicies(i)).last_abs;
         else
            cal_from = from_year;
   
         float log_p_uninv = 0;
         for(int t=cal_from+1;t<=2010;t++)
            log_p_uninv += -pow(alpha*lakes(indicies(i)).pp(t)+gamma_par,c_par);

          pred_p(m,i) = 1 - exp(log_p_uninv);
      }
      
   }
   val_sim_file.close();
   return(pred_p);
}





///////////////// Utilities ///////////////////

float average(_vbc_vec<float> x)
{
   int n=x.UBound();
   float tmp_sum=0;
   for(int i=1;i<=n;i++)
      tmp_sum += x(i);

   float nf=n;
   float avg= (1/nf)*tmp_sum;
   return(avg);
}
int which_min(int a,int b)
{
   if(a<=b)
       return a;
   else
       return b;
}
int which_max(int a,int b)
{
   if(a>=b)
       return a;
   else
       return b;
}
void write_t()
{
    for(int i=1;i<=n_lakes;i++)
        t_file<< t_vec(i) << "\t";

    t_file << "\n";
}
void write_par()
{
        for(int i=1;i<=n_chem_var+1;i++)
            par_file<< chem_pars(i) << "\t";

        par_file << d_par << "\t" << e_par;

    par_file << "\n";

    //for(int i=1;i<=n_lakes;i++)
    //    alpha_mcmc_file << calc_alpha(i) << "\t";

    //alpha_mcmc_file << "\n";
}
void write_traf_mat()
{
    ofstream traf;
    traf.open("output/traf_mat.gephi");
    traf << ";";
    for(int i=1;i<=n_lakes-1;i++)
        traf << "l" << i << ";";
    traf << "l" << n_lakes << "\n";

    for(int i=1;i<=n_lakes;i++)
    {
        traf << "l" << i << ";";
        for(int j=1;j<=n_lakes-1;j++)
            traf << traf_mat(i,j) << ";";
        traf << traf_mat(i,n_lakes) << "\n";
    }
    traf.close();
}
void write_inv_stat()
{
   ofstream inv_stat;
   inv_stat.open("output/inv_stat.dat");

   for(int i=1;i<=n_lakes;i++)
   {
      inv_stat << t_vec(i) << "\n";
   }

   inv_stat.close();
}
void write_pp()
{
    ofstream pp_file;
    pp_file.open("output/pp.dat");

    int lake_to_index;
    for(int t=from_year+1;t<=to_year;t++)
    {
        for(int j=1;j<=state(t).n_calc_pp;j++) //n_calc_pp=the number of uninvaded+newly invaded sites
        {
            lake_to_index=state(t).calc_pp_index(j);
            pp_file << t << "\t" << lake_to_index << "\t" << lakes(lake_to_index).pp(t)<<"\n";
        }
    }
   pp_file.close();
}

// Sample an integer vector with replacement. -- For generating bootstrap samples of lake indicies.
_vbc_vec<int> sample_w_replace(_vbc_vec<int> vec)
{
   _vbc_vec<int> s_vec(1,vec.UBound());
   int tmp_index;
   for(int i=1;i<=vec.UBound();i++)
   {
      tmp_index = (int) runif(1,vec.UBound()+1);
      s_vec(i) = vec(tmp_index);
   }

   return(s_vec);
}
// Sample an integer vector without replacement.
_vbc_vec<int> sample_wo_replace(_vbc_vec<int> vec,int n)
{
   int tmp_index,tmp;
   _vbc_vec<int> s_vec(1,n);
   //cout<< "In sample_wo_replace XX\n";
   for(int i=1;i<=n;i++)
   {
      tmp_index = (int) runif(1,vec.UBound() + 2 - i);
      tmp = vec(tmp_index);
      vec(tmp_index) = vec(vec.UBound() + 1 - i);
      vec(vec.UBound() + 1 - i) = tmp;
      s_vec(i)=tmp;
   }
   //cout<< "In sample_wo_replace YY\n";
   return(s_vec);
}

// Emulates wc -l (number of lines in a file)
int wc_l(string path_to_file)
{
   int number_of_lines=0;
   std::string line;
   std::ifstream myfile(path_to_file.c_str());
   while (std::getline(myfile, line))
        number_of_lines++;
   myfile.close();
   return(number_of_lines);
}

//Count the number of columns in a file
//**Only looks at first line, doesn't check sanity of entire file **
int wc_column(string path_to_file)
{
    int number_of_columns=0;
    string tmp;
    std::string line;
    std::ifstream myfile(path_to_file.c_str());
    while(!myfile.eof())
    {
        myfile >> tmp;
        number_of_columns++;
        if(myfile.peek() == '\n' ) break;
    }
    myfile.close();
    return(number_of_columns);
}

// given log(a) and log(b) return log(a + b)
float add_log(float la, float lb)
{
    return la + log(1 + exp(lb-la));
}

void normalize_vars()
{
    float sum(0);
    float mean;
    for(int j = 1;j<=n_chem_var;j++)
    {
        float sum(0), ss(0), sd;
        for(int i = 1; i<=n_lakes; i++)
            sum += lakes(i).chem(j);
        mean = sum/float(n_lakes);
        for(int i = 1; i<=n_lakes; i++)
            ss += pow(lakes(i).chem(j) - mean,2);
        sd = sqrt(ss/float(n_lakes));
        for(int i = 1; i<=n_lakes; i++)
            lakes(i).chem(j) = (lakes(i).chem(j) - mean) / sd;
    }
    for(int j = 1;j<=n_pdet_var;j++)
    {
        float sum(0), ss(0), sd;
        for(int i = 1; i<=n_lakes; i++)
            sum += lakes(i).pdet_var(j);
        mean = sum/float(n_lakes);
        for(int i = 1; i<=n_lakes; i++)
            ss += pow(lakes(i).pdet_var(j) - mean,2);
        sd = sqrt(ss/float(n_lakes));
        for(int i = 1; i<=n_lakes; i++)
            lakes(i).pdet_var(j) = (lakes(i).pdet_var(j) - mean) / sd;
    }
}



////////////////////////////////////////////////////////////////////////////////
///// UNCERTAINTY CALCULATIONS ON PROBABILITY OF ESTABLISHMENT /////////////////
////////////////////////////////////////////////////////////////////////////////
void save_beta_dists(_vbc_vec<float> * pr_inv_uncert)
{
    _vbc_vec<float> ab(1,2);
    for(int t=from_year;t<=to_year;t++)
    {   
        for(int i=1;i<=n_lakes;i++)
        {
            ab = fit_beta(pr_inv_uncert,i,t);
            cerr << ab(1) << " " << ab(2) << "\n";
        }
    } 
}

_vbc_vec<float> fit_beta(_vbc_vec<float> * pr_inv_uncert,int i, int t)
{
    int reps = (*pr_inv_uncert).UBound(1);
    _vbc_vec<float> dat1(1,1000);
    _vbc_vec<float> params1(1,2);
    _vbc_vec<float> MLE_params(1,2);
    params1(1) = 1;
    params1(2) = 1;
    for(int r=1;r<=1000;r++)
        dat1(r) = (*pr_inv_uncert)(r,i,t);

    simplex::clsSimplex<float> model_fit;
    model_fit.start(&dat1,&params1, &beta_lhood,2, 1e-10);
    model_fit.getParams(&MLE_params);
    return(MLE_params);
}
float beta_lhood(_vbc_vec<float> * pars, _vbc_vec<float> * dat)
{
    int n_dat = (*dat).UBound();
    float a=(*pars)(1);
    float b=(*pars)(2);
    if(a < 0 || b < 0)
        return(100000000);
    float ll(0);
    for(int i=1;i<=n_dat;i++)
    {
        ll += dbeta((*dat)(i),a,b,1);
    }   
    return(-ll);
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////




void print_lake(int i)
{
    cout << "\n\nlake " << i << " summary:\n";
    cout << "ID: " << lakes(i).id << "\n";
    cout << "Wj: " << lakes(i).Wj << "\n";
    cout << "Discovered: " << lakes(i).discovered << "\n";
    for(int t=from_year;t<=to_year;t++)
    {
        cout << t << ":" << lakes(i).pp(t) << " " << (1-exp(- pow(calc_alpha(i)  * lakes(i).pp(t) + gamma_par, c_par )))  << "\n";
    }
    cout << "\n";
    cout << "alpha: " << calc_alpha(i) << "\n";

 
    /*
		float Wj,x,y,alpha,Uj;
        bool invaded;
        bool val_invaded;
        int discovered;
        int last_abs;
        int status_2010;
        int id;
        _vbc_vec<float> chem;
	    _vbc_vec<float> Uij;
	    _vbc_vec<float> Qjt;
        _vbc_vec<float> pp;
        _vbc_vec<float> PInv;
        _vbc_vec<float> PInv_inst;
        _vbc_vec<float> pdet_var;
    */
}

// Parse params
// Read from a file defining the parameters to be used for fitting and simulating
void parse_par_seeds(string param_def_file,_vbc_vec<float> *parval)
{
    ifstream parseeds(param_def_file.c_str());
    int n_par = wc_l(param_def_file);
    (*parval).redim(1,n_par);
    parnames.redim(1,n_par);
    cout << "Seed params:\n"; 
    for(int i=1;i<=n_par;i++)
    {
        parseeds >> parnames(i);
        parseeds >> (*parval)(i);
        cout << "\t" << parnames(i) << "\t" << (*parval)(i) << "\n";        
    }
    parseeds.close();
    //Test
    //cout << glb_alpha << " :: " << d_par << " :: " << e_par << " :: " << c_par << "\n";
}


//Link par is called from parse_params
void link_pars(_vbc_vec<string> parname_str,_vbc_vec<float> parval)
{
    for(int i = 1; i<= parval.UBound(); i++)
    {
        if(parname_str(i) == "alpha")
            glb_alpha = exp(parval(i));//-log(1-(1/(1+exp(-parval(i)))));
        else if(parname_str(i) == "d")
        {
            one_time_traf = FALSE;
            traf_mat_precompute = TRUE;
            d_par = parval(i);
            glb_traf_array_ind = lookup_traf_mat_array_index();
        }
        else if(parname_str(i) == "e")
        {
            traf_mat_precompute = TRUE;
            e_par = parval(i);
            glb_traf_array_ind = lookup_traf_mat_array_index();
        }
        else if(parname_str(i) == "c")
            c_par = parval(i);
        else if(parname_str(i) == "gamma")
            gamma_par = parval(i);
        else if(parname_str(i) == "pdet")
        {
            pdet = parval(i);
            fit_pdet = TRUE;
            hetero_pdet = FALSE;
        }
        else if(parname_str(i) == "pdet_var0")
        {
            est_pdet = 0;
            pdet_pars(1) = parval(i);
            fit_pdet = TRUE;
            hetero_pdet = TRUE;
        }
        else if(parname_str(i) == "pdet_var1")
        {
            est_pdet = 1;
            pdet_pars(2) = parval(i);
            fit_pdet = TRUE;
            hetero_pdet = TRUE;
        }
        else if(parname_str(i) == "pdet_var2")
        {
            est_pdet = 2;
            pdet_pars(3) = parval(i);
            fit_pdet = TRUE;
            hetero_pdet = TRUE;
        }
        else if(parname_str(i) == "pdet_var3")
        {
            est_pdet = 3;
            pdet_pars(4) = parval(i);
            fit_pdet = TRUE;
            hetero_pdet = TRUE;
        }
        else if(parname_str(i) == "pdet_var4")
        {
            est_pdet = 4;
            pdet_pars(5) = parval(i);
            fit_pdet = TRUE;
            hetero_pdet = TRUE;
        }
        else if(parname_str(i) == "env0") //intercept
        {
            chem_pars(1) = parval(i);
            env = TRUE;
        }
        else if(parname_str(i) == "env1"){
            est_env = 1;
            chem_pars(2) = parval(i);
        }
        else if(parname_str(i) == "env2"){
            est_env = 2;
            chem_pars(3) = parval(i);
        }
        else if(parname_str(i) == "env3"){
            est_env = 3;
            chem_pars(4) = parval(i);
        }
        else if(parname_str(i) == "env4"){
            est_env = 4;
            chem_pars(5) = parval(i);
        }
        else if(parname_str(i) == "env5"){
            est_env = 5;
            chem_pars(6) = parval(i);
        }
        else if(parname_str(i) == "env6"){
            est_env = 6;
            chem_pars(7) = parval(i);
        }
        else if(parname_str(i) == "env7"){
            est_env = 7;
            chem_pars(8) = parval(i);
        }
        else if(parname_str(i) == "env8"){
            est_env = 8;
            chem_pars(9) = parval(i);
        }
        else if(parname_str(i) == "env9"){
            est_env = 9;
            chem_pars(10) = parval(i);
        }
        else if(parname_str(i) == "env10"){
            est_env = 10;
            chem_pars(11) = parval(i);
        }
        else if(parname_str(i) == "env11"){
            est_env = 11;
            chem_pars(12) = parval(i);
        }
        else if(parname_str(i) == "env12"){
            est_env = 12;
            chem_pars(13) = parval(i);
        }
        else if(parname_str(i) == "env13"){
            est_env = 13;
            chem_pars(14) = parval(i);
        }
        else if(parname_str(i) == "env14"){
            est_env = 14;
            chem_pars(15) = parval(i);
        }
        else if(parname_str(i) == "env15"){
            est_env = 15;
            chem_pars(16) = parval(i);
        }
        else if(parname_str(i) == "env16"){
            est_env = 16;
            chem_pars(17) = parval(i);
        }
        else if(parname_str(i) == "env17"){
            est_env = 17;
            chem_pars(18) = parval(i);
        }
        else if(parname_str(i) == "env18"){
            est_env = 18;
            chem_pars(19) = parval(i);
        }
        else if(parname_str(i) == "env19"){
            est_env = 19;
            chem_pars(20) = parval(i);
        }
        else if(parname_str(i) == "env20"){
            est_env = 20;
            chem_pars(21) = parval(i);
        }

        // 2nd order terms ** These terms must be added AFTER all first order terms
        // and must be added in order //
        else if(parname_str(i) == "env1_2"){
            alpha_order = 2;
            chem_pars(est_env+2) = parval(i);
        }
        else if(parname_str(i) == "env2_2"){
            alpha_order = 2;
            chem_pars(est_env+3) = parval(i);
        }
        else if(parname_str(i) == "env3_2"){
            alpha_order = 2;
            chem_pars(est_env+4) = parval(i);
        }
        else if(parname_str(i) == "env4_2"){
            alpha_order = 2;
            chem_pars(est_env+5) = parval(i);
        }
        else if(parname_str(i) == "env5_2"){
            alpha_order = 2;
            chem_pars(est_env+6) = parval(i);
        }
        else if(parname_str(i) == "env6_2"){
            alpha_order = 2;
            chem_pars(est_env+7) = parval(i);
        }
        else if(parname_str(i) == "env7_2"){
            alpha_order = 2;
            chem_pars(est_env+8) = parval(i);
        }
        else if(parname_str(i) == "env8_2"){
            alpha_order = 2;
            chem_pars(est_env+9) = parval(i);
        }
        else if(parname_str(i) == "env9_2"){
            alpha_order = 2;
            chem_pars(est_env+10) = parval(i);
        }
        else if(parname_str(i) == "env10_2"){
            alpha_order = 2;
            chem_pars(est_env+11) = parval(i);
        }

        else if(parname_str(i) == "dummy")
            float tmp = parval(i);
        else
        {
            cout << "Malformed parameter name\n"; 
            exit(1);
        }
    }
}











////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
