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
#include"man_grav.h"

// Read in sim inits //
int main(int argc,char *argv[])
{
    // decode arguments
    args(argc,argv);
    set_seed(rngseed,54321);
    
    if(output_folder != "")
        int status = mkdir(output_folder.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);


    // Confusingly, the output_folder is actually the input folder (scientific computing!!!)
    // so we'll save that here so that we can read in the mcmc file.
    string mcmc_in_folder = output_folder;

    // Then we'll add the management output folder as a sub-folder of the input
    output_folder = output_folder + "/management";
    int status2 = mkdir(output_folder.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    impose_management = TRUE;
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

    // TEST //
    //_vbc_vec<float> pr_inv_2013a(1,n_lakes);
    //for(int i=1;i<=n_lakes;i++)
    //    pr_inv_2013a(i) = runif(0,1);
    //set_top_n_source_lakes(500, 1, 1, pr_inv_2013a);
    // -- //



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
      stringstream out;
      out << r;
      glb_r = r; 
      string s = out.str();

      string t_vec_filename = output_folder + "/tvec" + s;
      ofstream t_vec_file(t_vec_filename.c_str());
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


    //////////////// MCMC ////////////////////////////
    if(run_type==2){
        //force_invasions = TRUE;
        string mcmc_file;
        mcmc_file = mcmc_in_folder + "/sim_lib.mcmc1";
        cout << mcmc_file << "\n";

        int n_pars = params1.UBound();

        //////////////// sim inv given posterior ////////////
        //compensate_force = FALSE;
        //only_known = TRUE;

        ifstream posterior_file(mcmc_file.c_str());
        int n_mcmc = wc_l(mcmc_file);
        int burn_in = 1000;
        float tmp;
        to_year = 2023;
        _vbc_vec<float> pr_inv(1,n_lakes,from_year,to_year);
        for(int i=1;i<=n_lakes;i++)
        {
            for(int t=from_year;t<=to_year;t++)
            {
                pr_inv(i,t) = 0;
            }
        }
        //force_invasions = TRUE;

        // First run to find pr_inv for 2013 (no management)
        impose_management = FALSE;
        //n_mcmc = 2000;
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
                cout << float(state(to_year).n_inv)/float(n_lakes) << "\t" <<  float(state(2012).n_inv)/float(n_lakes) << "\n";;
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
                }
                // ----- //
            }
        }
        posterior_file.close();


        // Now apply management at top n outflow sites //
        impose_management = TRUE;

        _vbc_vec<float> pr_inv_2013(1,n_lakes);
        for(int i=1;i<=n_lakes;i++)
        {
            pr_inv_2013(i) = pr_inv(i,2013);
            //cout << pr_inv_2013(i) << " :: ";
        }


        // reset pr_inv //        
        for(int i=1;i<=n_lakes;i++)
        {
            for(int t=from_year;t<=to_year;t++)
            {
                pr_inv(i,t) = 0;
            }   
        }
        // ----------- //

        ////////////////////////////////////////////////////////////////////////
        // Management Scenarios //
        ////////////////////////////////////////////////////////////////////////
        int topn;
        _vbc_vec<int> topn_vec(1,8);
        _vbc_vec<int> topn_index; // index of the top n lakes selected for management
        topn_vec(1) = 1;
        topn_vec(2) = 2;
        topn_vec(3) = 4;
        topn_vec(4) = 8;
        topn_vec(5) = 16;
        topn_vec(6) = 32;
        topn_vec(7) = 64;
        topn_vec(8) = 128;
        float theta;
        for(int thetas=1;thetas<=2;thetas++)
        {
            if(thetas==1)
                theta = 0.35; // at full cost to boaters

            if(thetas==2)
                theta = 0.97; // at no cost to boaters

            for(int topn_count=1;topn_count<=8;topn_count++)
            {
                topn = topn_vec(topn_count);
                stringstream outname;
                outname << "_" << topn << "_" << theta;
                string nameadd = outname.str();

                topn_index = set_top_n_source_lakes(topn, theta, 0.01, pr_inv_2013); //(ntop,theta,phi,pr_inv)
                calc_traf_mat_array_man();
                // ------------------------------------------ //

                string pred_file_string;
                pred_file_string = output_folder + "/pred_prob" + nameadd + ".csv";
                ofstream pred_file(pred_file_string.c_str());

                string traject_file_string;
                traject_file_string = output_folder + "/traject" + nameadd + ".csv";
                ofstream traject_file(traject_file_string.c_str());

                string topn_file_string;
                topn_file_string = output_folder + "/topn" + nameadd + ".csv";
                ofstream topn_file(topn_file_string.c_str());
                for(int n=1;n<=topn;n++)
                    topn_file << topn_index(n) << "\n";


                // This time management kicks in in 2013 //
                posterior_file.open(mcmc_file.c_str());
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
                        for(int t=2013;t<=to_year;t++)
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
                    }
                }
                // Mean probabilities over time //
                for(int t=2013;t<=to_year;t++)
                {
                    for(int i=1;i<=n_lakes-1;i++)
                    {
                        pred_file << pr_inv(i,t) << ",";
                        pr_inv(i,t) = 0; //reset //
                    }
                    pred_file << pr_inv(n_lakes,t) << "\n";
                    pr_inv(n_lakes,t) = 0; //reset //
                }

                pred_file.close();
                traject_file.close();
                posterior_file.close();
                topn_file.close();
                //compensate_force = TRUE;
                // //
            }
        }
    }
    return 0;
}



