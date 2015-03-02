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

// Read in sim inits //
int main(int argc,char *argv[])
{
    // decode arguments
    args(argc,argv);
    set_seed(rngseed,54321);
    if(output_folder != "")
        int status = mkdir(output_folder.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    inits();
    to_year = 2023;
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
    to_year = 2013;
    int n_inv = 0;
    int n_seed = 0;
    for(int i=1;i<=n_lakes;i++)
    {      
        if(t_vec(i) <= to_year){n_inv++;}
        if(t_vec(i) == from_year){n_seed++;}
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
            2000,
            100,
            1,
            mcmc_file.c_str(),
            true,
            true,
            true,
            200,
            4);

        //////////////// sim inv given posterior ////////////
        //compensate_force = FALSE;
        //only_known = TRUE;
        string pred_file_string;
        pred_file_string = output_folder + "/pred_prob" + s + ".csv";
        ofstream pred_file(pred_file_string.c_str());

        string pred_file_string_uncert2013;
        pred_file_string_uncert2013 = output_folder + "/pred_prob_uncert2013-" + s + ".csv";
        ofstream pred_file_uncert2013(pred_file_string_uncert2013.c_str());

        string pred_file_string_uncert2023;
        pred_file_string_uncert2023 = output_folder + "/pred_prob_uncert2023-" + s + ".csv";
        ofstream pred_file_uncert2023(pred_file_string_uncert2023.c_str());

        string traject_file_string;
        traject_file_string = output_folder + "/traject" + s + ".csv";
        ofstream traject_file(traject_file_string.c_str());


        string pr_inv_not_det_file_string;
        pr_inv_not_det_file_string = output_folder + "/pr_inv_not_det" + s + ".csv";
        ofstream pr_inv_not_det_file(pr_inv_not_det_file_string.c_str());


        ifstream posterior_file(mcmc_file.c_str());
        int n_mcmc = wc_l(mcmc_file);
        int burn_in = 0; //1000
        float tmp;
        to_year = 2023;
        _vbc_vec<float> pr_inv(1,n_lakes,from_year,to_year);
        _vbc_vec<float> pr_inv_uncert(1,n_mcmc-burn_in,1,n_lakes,1,2);
        for(int i=1;i<=n_lakes;i++)
        {
            for(int t=from_year;t<=to_year;t++)
            {
                pr_inv(i,t) = 0;
            }
            for(int t=1;t<=2;t++)
            {
                for(int rep=1;rep<=n_mcmc-burn_in;rep++)
                {
                    pr_inv_uncert(rep,i,t) = 0;
                }
            }
        }
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
                cout << float(state(to_year).n_inv)/float(n_lakes) << "\n";
                // predicted prob across time //
                for(int t=from_year;t<=to_year;t++)
                {   
                    int n_inv(0);
                    for(int i=1;i<=n_lakes;i++)
                    {
                        if(t_vec(i) <= t) 
                        {
                            pr_inv(i,t) += 1.0/float(n_mcmc-burn_in);
                            n_inv++;
                        }
                    }
                    if(t < to_year)
                        traject_file << float(n_inv)/float(n_lakes) << ",";
                    else
                        traject_file << float(n_inv)/float(n_lakes) << "\n";
                }
                // ----- //

                // UNCERTAINTY IN PROBABILITY OF ESTABLISHMENT //
                /*
                int t;
                for(int rep=1;rep<=1000;rep++)
                {
                    sim_spread();
                    // predicted prob across time for each sample from the posterior //
                    for(int ti=1;ti<=2;ti++)
                    {   
                        if(ti==1){t=2013;}
                        if(ti==2){t=2023;}
                        for(int i=1;i<=n_lakes;i++)
                        {
                            if(t_vec(i) <= t)
                                pr_inv_uncert(mc-burn_in,i,ti) += 1/(1000);
                        }
                    }
                }
                */
                // Prob invaded but not detected by 2013 //
                for(int i=1;i<=n_lakes-1;i++)
                    pr_inv_not_det_file << pr_inv_not_det(i,2013) << ",";
                pr_inv_not_det_file << pr_inv_not_det(n_lakes,2013) << "\n";
            }
        }
        // Mean probabilities over time //
        for(int t=from_year;t<=to_year;t++)
        {
            for(int i=1;i<=n_lakes-1;i++)
            {
                pred_file << pr_inv(i,t) << ",";
            }
            pred_file << pr_inv(n_lakes,t) << "\n";
        }

        // Uncertainty at 2013 and 2023 //
        for(int mc=1;mc<=n_mcmc-burn_in;mc++)
        {
            for(int i=1;i<=n_lakes-1;i++)
                pred_file_uncert2013 << pr_inv_uncert(mc,i,1) << ",";
            pred_file_uncert2013 << pr_inv_uncert(mc,n_lakes,1) << "\n";

            for(int i=1;i<=n_lakes-1;i++)
                pred_file_uncert2023 << pr_inv_uncert(mc,i,2) << ",";
            pred_file_uncert2023 << pr_inv_uncert(mc,n_lakes,2) << "\n";
        }
        pr_inv_not_det_file.close();
        pred_file_uncert2013.close();
        pred_file_uncert2023.close();
        pred_file.close();
        traject_file.close();
        posterior_file.close();
        //compensate_force = TRUE;
        // //
    }
    return 0;
}



