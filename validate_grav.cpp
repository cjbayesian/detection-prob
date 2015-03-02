//////////////////////////////////////////////////////////
// Simulates and Tests a Bayesian implimentation of 
// gravity mediated spread of invasive species
// using a presence only observation model
// 
// 	Corey Chivers, 2014
// 	McGill University
//	Department of Biology
//	corey.chivers@mail.mcgill.ca
//
// See Leung et al. 2006 for underlying 
// gravity and likelihood equation defs
// 
//////////////////////////////////////////////////////////

#include"grav.h"
//#include"sim_grav.h"

//returns the probability that site <lake_index> will be detected by <to_year>
// validation compares this probability to actuall dectections occuring after the training period.
float pr_det(int lake_index,int val_year)
{
    float tmp_lh;
    float pr_undet=0;
    float alpha = calc_alpha(lake_index);
    float pdet_ = calc_pdet(lake_index);
    for(int t_inv=from_year+1; t_inv<=to_year+1; t_inv++)
    {
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
        pr_undet += tmp_lh;
    }
    // Add in all of the probability of detections happening before val_year
    // since we know this did not happen and therefore what we  
    // want is P(detection in year > val_year | no dectection < val_year)
    for(int disc=from_year+1;disc<=val_year;disc++)
    {
        for(int t_inv=from_year+1;t_inv<=disc;t_inv++)
        {
            tmp_lh = 1; // The summed probability across states //
            for(int t=from_year+1;t<=disc;t++)
            {
                if(t < t_inv)
                    tmp_lh *= exp( -pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par) );
                if(t == t_inv)
                    tmp_lh *= 1 - exp( -pow(alpha*lakes(lake_index).pp(t)+gamma_par,c_par) ); 
                if(t >= t_inv && t != disc)
                    tmp_lh *= 1 - pdet_;
            }
            tmp_lh *= pdet_;
            pr_undet += tmp_lh;
        }        
    }
    return(1-pr_undet);
}


int main(int argc,char *argv[])
{
    // decode arguments
    args(argc,argv);
    set_seed(rngseed,54321);
    int status(0);
    if(output_folder != "")
        status = mkdir(output_folder.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if(status != 0)
    {
        cerr << "Could not create output folder (path to folder root must exist)\n";
        exit(0);
    }
    inits();
    to_year = 2013;
    int val_year = 2005;
    int r = 1;
    // Read in Data // 
    cout << "Reading Spatial Data\n";
    if(sites_file_name == "../2010_bytho_data/lakes_processed_normalized.csv")
        read_data();
    else
        read_po_data();


    // Set initial Params //
    cout << "Initiallizing...\n";
    e_par = 0.5;
    d_par = 2;
    gamma_par = 0;
    n_sampled = 0;



    clear_traf_mat();
    int n_pars; 
    int n_disc;
    _vbc_vec<float> params1;

    string pars_filename;
    pars_filename = output_folder + "/pars.csv";
    ofstream parnames_file(pars_filename.c_str());
    parse_par_seeds(parseeds_file,&params1);
    link_pars(parnames,params1);
    n_pars = params1.UBound();

    for(int p=1;p<=n_pars-1;p++)
        parnames_file << parnames(p) <<",";
    parnames_file << parnames(n_pars) << "\n";
    parnames_file.flush();
    parnames_file.close();

    init_state();
    init_t(); 
    calc_state();

    if(traf_mat_precompute)
    {
        calc_traf_mat_array();
        //read_traf_mat_array();
    }
    else {
        calc_traf();
        calc_Uj();
        calc_traf_mat();
        //write_traf_mat();
    }
    parse_par_seeds(parseeds_file,&params1);
    to_year = val_year; //train on discoveries up to this year.
    _vbc_vec<int> all_discoveries(1,n_lakes);    
    int n_inv = 0;
    int n_seed = 0;
    for(int i=1;i<=n_lakes;i++)
    {      
        if(t_vec(i) <= to_year){n_inv++;}
        if(t_vec(i) == from_year){n_seed++;}
        all_discoveries(i) = lakes(i).discovered; //save all discoveries (training and test period)
        if(lakes(i).discovered != 0 && lakes(i).discovered > to_year) //remove any discoveries after the training period
        {
            lakes(i).discovered = 0;
        }
    }
    
    cout << n_seed << " seed lakes\n";
    gen_prop = float(n_inv)/float(n_lakes); //proportion of sites detected

    n_disc = 0;
    for(int i=1;i<=n_lakes;i++)
        if(lakes(i).discovered != 0 && lakes(i).discovered != from_year) {n_disc++;}
    cout << n_disc << " lakes with discoveries.\n";

    which_sampled_or_valid();
    //////////////////////////////////////////////////
      ofstream t_vec_file("tvec.csv");
      for(int i=1;i<=n_lakes-1;i++)
         t_vec_file << t_vec(i) << ",";
      t_vec_file << t_vec(n_lakes) << "\n";
      
      for(int i=1;i<=n_lakes-1;i++)
         t_vec_file << lakes(i).discovered << ",";
      t_vec_file << lakes(n_lakes).discovered << "\n";

      for(int i=1;i<=n_lakes-1;i++)
         t_vec_file << lakes(i).last_abs << ",";
      t_vec_file << lakes(n_lakes).last_abs << "\n";

      for(int i=1;i<=n_lakes-1;i++)
         t_vec_file << lakes(i).Wj << ",";
      t_vec_file << lakes(n_lakes).Wj << "\n";

      for(int i=1;i<=n_lakes-1;i++)
         t_vec_file << lakes(i).x << ",";
      t_vec_file << lakes(n_lakes).x << "\n";

      for(int i=1;i<=n_lakes-1;i++)
         t_vec_file << lakes(i).y << ",";
      t_vec_file << lakes(n_lakes).y << "\n";

      for(int i=1;i<=n_lakes-1;i++)
         t_vec_file << from_year << ",";
      t_vec_file << 2023 << "\n";
      //////////////////////////////////////////////////
      t_vec_file.close();
      stringstream out;
      out << r;
      glb_r = r; 
      string s = out.str();

      string cp_cmd("cp tvec.csv " + output_folder + "/tvec");
      cp_cmd = cp_cmd + s;
      system(cp_cmd.c_str());

    //////////////// MCMC ////////////////////////////
    if(run_type==2){
        //force_invasions = TRUE;
        string mcmc_file;
        mcmc_file = output_folder + "/sim_lib.mcmc" + s;
        cout << mcmc_file << "\n";

        _vbc_vec<float> prop_width;

        int n_pars = params1.UBound();
        prop_width.redim(1,n_pars);

        for(int p=1;p<=n_pars;p++){
            prop_width(p)=0.005;
        }

        _vbc_vec<float> prop_sigma;
        prop_sigma = diag(prop_width);

        mcmcMD::run_mcmc(params1, 
            prop_sigma,
            &likelihood_wrapperMCMC_MD,
            &prior_MD,
            &restrict_MCMC_MD, 
            30000,
            500,
            4,
            mcmc_file.c_str(),
            true,
            true,
            true,
            500,
            5);

        //////////////// sim inv given posterior ////////////
        //compensate_force = FALSE;
        //only_known = TRUE;
        string pred_file_string;
        pred_file_string = output_folder + "/pred_det_prob" + s + ".csv";
        ofstream pred_file(pred_file_string.c_str());
        for(int i=1;i<=n_lakes-1;i++)
        {
             lakes(i).discovered = all_discoveries(i); //add test discoveries
             pred_file << lakes(i).discovered << ",";
        }
        lakes(n_lakes).discovered = all_discoveries(n_lakes); //add test discoveries
        pred_file << lakes(n_lakes).discovered << "\n";

        string pred_alpha_file_string;
        pred_alpha_file_string = output_folder + "/pred_alpha" + s + ".csv";
        ofstream pred_alpha_file(pred_alpha_file_string.c_str());
        

        ifstream posterior_file(mcmc_file.c_str());
        int n_mcmc = wc_l(mcmc_file);
        int burn_in = 1000;
        float tmp;
        to_year = 2013;

        //force_invasions = TRUE;
        for(int mc=1;mc<=n_mcmc;mc++)
        {  
            posterior_file >> tmp; //index
            for(int pa=1;pa<=n_pars;pa++)
                posterior_file >> params1(pa);
            posterior_file >> tmp; //log-likelihood
            if(mc > burn_in) //discard burn_in
            {
                link_pars(parnames,params1);
                sim_spread();
                for(int i=1;i<=n_lakes-1;i++){
                    pred_file << pr_det(i,val_year) << ",";
                    pred_alpha_file << calc_alpha(i) << ",";
                }
                pred_file << pr_det(n_lakes,val_year) << "\n";
                pred_alpha_file << calc_alpha(n_lakes) << "\n";
                cout << float(state(to_year).n_inv)/float(n_lakes) << "\n";
            }
        }
        pred_file.close();
        posterior_file.close();
        pred_alpha_file.close();
        // //
    }
    return 0;
}



