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
#include"sim_grav.h"

// Read in sim inits //
int main(int argc,char *argv[])
{
   // decode arguments
   args(argc,argv);
   set_seed(rngseed,54321);
   if(output_folder != "")
      int status = mkdir(output_folder.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

   // Sim lakes and sources //
   // <<< Instead just use 2EB or Ontario lakes
   inits();
   to_year = 2023;
   sim_inits();
   //system("mkdir test_folder");
   // Read in Data // 
   cout << "Reading Spatial Data\n";
   if(sites_file_name == "../2010_bytho_data/lakes_processed_normalized.csv")
        read_data();
   else
        read_po_data();


   // Set initial Params //
   cout << "Initiallizing...\n";
   // Initiallize (seed) the largest lake //
   // <<< Instead just use actual seed lakes
   // <<< But we will need to clear the current records of
   //     last_abs, last_obs, last_obs_uninv, discovered

   //for faster sims:
   //n_lakes = 500;
   n_sources = 100;
   clear_traf_mat();
   

   //one_time_traf = TRUE;
   traf_mat_precompute = TRUE;
   calc_traf_mat_array();
   //read_traf_mat_array();
  
   sim = TRUE;
   int n_pars; 
   int n_disc;
   _vbc_vec<float> true_params;
   _vbc_vec<float> params1;
   _vbc_vec<float> dat1;
   _vbc_vec<float> MLE_params;
   string mle_filename;
   mle_filename = output_folder + "/mles.csv";
   ofstream mle_pars(mle_filename.c_str());
   //ofstream mean_field_file("output/mean_field.csv");
   parse_par_seeds(parseeds_file,&params1);
   n_pars = params1.UBound();
   for(int p=1;p<=n_pars-1;p++)
       mle_pars << parnames(p) <<",";
   mle_pars << parnames(n_pars) << "\n";
   mle_pars.flush();
   //for(int i=1;i<=n_lakes;i++)
   //    lakes(i).Wj = log(lakes(i).Wj);
   for(int r=1;r<=100;r++)
   {
      n_sampled = 0;
      // Drop situations where very few lakes are observed
      //only_known = FALSE;
      while(n_sampled < 15)
      {
         to_year = 2023;
         for(int i=1;i<=n_lakes;i++)
         {
            //clear all but the seed lake(s)
            lakes(i).last_abs = 0;
            lakes(i).status_2010 = 0;
            if(lakes(i).discovered != from_year)
            {
               lakes(i).discovered = 0;
               lakes(i).invaded = 0;
            }
         }
         init_state();
         init_t(); 
         calc_state();
         generate_params();

         // Simulate the invasion //
         //calc_traf();
         //calc_Uj();
         clear_t();
         //if(one_time_traf)
         //   calc_traf_mat();
         //write_traf_mat();

         sim_spread();
         to_year = 2013;
         int n_inv = 0;
         int n_seed = 0;
         for(int i=1;i<=n_lakes;i++)
         {      
            if(t_vec(i) <= to_year){n_inv++;}
            if(t_vec(i) == from_year){n_seed++;}
         }         
         cout << n_inv << " invaded of " << n_lakes << " invaded\n";
         cout << n_seed << " seed lakes\n";
         if(n_inv == n_seed || n_inv == n_lakes) {continue;}
         gen_prop = float(n_inv)/float(n_lakes);

         true_params = link_gen_pars(parnames);
         glb_true_params = true_params;
         // Simulate the detection process //
         if(fit_pdet)
            detect();
         else
            sample_PA();


         n_disc = 0;
         for(int i=1;i<=n_lakes;i++)
            if(lakes(i).discovered != 0 && lakes(i).discovered != from_year) {n_disc++;}
         cout << n_disc << " lakes with discoveries.\n";
         if(n_disc == 0) {continue;}

         which_sampled_or_valid();
      }

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
      //////////////////////////////////////////////////
      t_vec_file.close();

      
      //////////////// MLE /////////////////////////////
      if(run_type==1){
         n_pars = params1.UBound();
         dat1.redim(1,n_pars);
         MLE_params = params1;

         simplex::clsSimplex<float> model_fit;
         model_fit.start(&dat1,&params1, &MLE_l_hood,n_pars, 1e-10);
         model_fit.getParams(&MLE_params);

         cout << "\n";
         for(int p=1;p<=n_pars;p++)
            cout << setw(15) << parnames(p) <<"\t";
         cout << "\n";
         for(int p=1;p<=n_pars;p++)
            cout << setw(15) << MLE_params(p) <<"\t";
         cout << "\n";
         for(int p=1;p<=n_pars;p++)
            cout << setw(15) << true_params(p) <<"\t";
         cout << "\n";
         
         // Write MLE and true params to file //
         for(int p=1;p<=n_pars;p++)
            mle_pars << MLE_params(p) <<",";
         for(int p=1;p<=n_pars-1;p++)
            mle_pars << true_params(p) <<",";
         mle_pars << true_params(n_pars) << "\n";
         mle_pars.flush();
      }


    //////////////// MCMC ////////////////////////////
    if(run_type==2){
        //only_known = TRUE;
        //force_invasions = TRUE;
        string mcmc_file;
        mcmc_file = output_folder + "/sim_lib.mcmc" + s;
        gen_pdet = pdet;
        cout << mcmc_file << "\n";

        _vbc_vec<float> prop_width;

        int n_pars = params1.UBound();
        prop_width.redim(1,n_pars);

        for(int p=1;p<=n_pars-1;p++)
            mle_pars << true_params(p) <<",";
        mle_pars << true_params(n_pars) << "\n";
        mle_pars.flush();

        params1 = true_params; //seed with generating values (give a leg-up)
float tt;        
likelihood_wrapperMCMC_MD(&params1,&tt,r);

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
            2,
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

        string traject_file_string;
        traject_file_string = output_folder + "/traject" + s + ".csv";
        ofstream traject_file(traject_file_string.c_str());

        ifstream posterior_file(mcmc_file.c_str());
        int n_mcmc = wc_l(mcmc_file);
        float tmp;
        to_year = 2023;
        _vbc_vec<float> pr_inv(1,n_lakes,from_year,to_year);
        for(int i=1;i<=n_lakes;i++)
        {
            for(int t=from_year;t<=to_year;t++)
                pr_inv(i,t) = 0;
        }
        //force_invasions = TRUE;
        int nn(0);
        for(int mc=1;mc<=n_mcmc;mc++)
        {  
            posterior_file >> tmp; //index
            for(int pa=1;pa<=n_pars;pa++)
                posterior_file >> params1(pa);

            posterior_file >> tmp; //log-likelihood
            if(mc > 400) //discard burn_in
            {
                link_pars(parnames,params1);
                sim_spread();
                pr_inv = prob_inv();
                //pr_inv_file << glb_r << "," << gen_prop << "," << float(state(to_year).n_inv)/float(n_lakes) << "\n";
                cout << glb_r << "," << gen_prop << "," << float(state(to_year).n_inv)/float(n_lakes) << "\n";
                // predicted prob across time //
                for(int t=from_year+1;t<=to_year;t++)
                {   
                    //int n_inv(0);
                    if(t == 2013 || t == 2018 || t == 2023)
                    {
                        pred_file << t << ",";
                        for(int i=1;i<=n_lakes-1;i++)
                        {
                            pred_file << pr_inv(i,t) << ",";
                            /*
                            if(t_vec(i) <= t) 
                            {
                                pr_inv(i,t) = pr_inv(i,t) + 1;
                                n_inv++;
                            }
                            */
                        }
                        pred_file << pr_inv(n_lakes,t) << "\n";
                    }
                    
                    //if(t < to_year)
                    //    traject_file << float(n_inv)/float(n_lakes) << ",";
                    //else
                    //    traject_file << float(n_inv)/float(n_lakes) << "\n";
                }
                nn++;
                // ----- //
            }
        }
        /*
        for(int t=from_year;t<=to_year;t++)
        {
            for(int i=1;i<=n_lakes-1;i++)
            {
                pred_file << pr_inv(i,t)/nn << ",";
            }
            pred_file << pr_inv(n_lakes,t)/nn << "\n";
        } 
        */
        pred_file.close();
        traject_file.close();
        posterior_file.close();
        //compensate_force = TRUE;
        // //
    }
    if(run_type == 3)
    {
        params1 = true_params; //seed with generating values (give a leg-up)
        float tt;
        for(int i=1;i<=40;i++)
            likelihood_wrapperMCMC_MD(&params1,&tt,r);
    }

   }
   //mean_field_file.close();
   mle_pars.close();
   pr_inv_file.close(); 
   return 0;
}


