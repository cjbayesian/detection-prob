#ifndef SIM_GRAV_H
#define SIM_GRAV_H

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  SIMS  ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Globals //
   float glb_alpha_min, glb_alpha_max, c_par_min, c_par_max;
   float d_par_min, d_par_max, e_par_min, e_par_max, pdet_min, pdet_max;
   float env0_min, env0_max, env_min, env_max;
   float pdet0_max, pdet0_min, pdet_c_min, pdet_c_max;

void sim_inits();
void generate_params();
_vbc_vec<float> link_gen_pars(_vbc_vec<string>, _vbc_vec<float>);
_vbc_vec<float> prob_inv();

#endif;


// Read in sim inits //
void sim_inits()
{
   string tmp;
   ifstream init_file("inits_sim.ini");

   init_file >> tmp;
   init_file >> glb_alpha_min;
   init_file >> glb_alpha_max;
   init_file >> tmp;
   init_file >> c_par_min;
   init_file >> c_par_max;
   init_file >> tmp;
   init_file >> d_par_min;
   init_file >> d_par_max;
   init_file >> tmp;
   init_file >> e_par_min;
   init_file >> e_par_max;
   init_file >> tmp;
   init_file >> pdet_min;
   init_file >> pdet_max;
   init_file >> tmp;
   init_file >> env0_min;
   init_file >> env0_max;
   init_file >> tmp;
   init_file >> env_min;
   init_file >> env_max;   
   init_file >> tmp;
   init_file >> pdet0_min;
   init_file >> pdet0_max; 
   init_file >> tmp;
   init_file >> pdet_c_min; //coeff on pdet, not pdet itself.
   init_file >> pdet_c_max; 

   init_file.close();
}

void generate_params()
{
   glb_alpha = runif(glb_alpha_min,glb_alpha_max);
   c_par = runif(c_par_min,c_par_max);
   d_par = runif(d_par_min,d_par_max);
   e_par = runif(e_par_min,e_par_max);
   pdet = runif(pdet_min,pdet_max);
   chem_pars(1) = runif(env0_min,env0_max);
   for(int ch=2;ch<=est_env+1;ch++)
      chem_pars(ch) = runif(env_min,env_max);
   pdet_pars(1) = runif(pdet0_min,pdet0_max);
   for(int ch=2;ch<=est_pdet+1;ch++)
      pdet_pars(ch) = runif(pdet_c_min,pdet_c_max);

   glb_traf_array_ind = lookup_traf_mat_array_index();
    cout << "TRAF: d=" << d_par << " e=" << e_par << " ind:" << glb_traf_array_ind << "\n";
   //cout << glb_alpha << "\t" << c_par << "\t" << d_par << "\t" << e_par << "\t" << pdet << "\n";
}

void detect()
{
   // Which sites to monitor
   int n_monitored = n_lakes;
   _vbc_vec<int> lake_inds(1,n_lakes);
   for(int i=1;i<=n_lakes;i++) {lake_inds(i)=i;}
   //lake_inds = sample_wo_replace(lake_inds,n_monitored);
   bool assume_pdet_1 = FALSE;
   float pdet_;
   // Monitor until detected or to_year
   int lake_index;
   for(int i=1;i<=n_monitored;i++)
   {
      lake_index = lake_inds(i);
      pdet_ = calc_pdet(lake_index);
      for(int t=t_vec(lake_index);t<=to_year;t++)
      {
         if(runif(0,1) < pdet_ && lakes(lake_index).discovered == 0)
         {
            //cout << t_vec(lake_index) << "\t" << t << "\n"; 
            lakes(lake_index).discovered = t;
            //if(assume_pdet_1)
            //   lakes(lake_index).last_abs = t-1;
            break;
         }
      }

   }
   //for(int i=1;i<=n_lakes;i++)
   //   cout << t_vec(i) << " :: detected :: " << lakes(i).discovered << "\n";
}

void sample_PA()
{
   // Which sites to monitor
   int n_monitored = n_lakes;
   _vbc_vec<int> lake_inds(1,n_lakes);
   for(int i=1;i<=n_lakes;i++) {lake_inds(i)=i;}
   lake_inds = sample_wo_replace(lake_inds,n_monitored);
   
   int n_sample_years = 4;
   _vbc_vec<int> sample_years(1,n_sample_years);
   sample_years(1) = 1995;
   sample_years(2) = 1999;
   sample_years(3) = 2003;
   sample_years(4) = 2005;
   float pr_sample = 0.2;

   for(int t=1;t<=n_sample_years;t++)
   {
      for(int i=1;i<=n_monitored;i++)
      {
         if(t_vec(lake_inds(i)) <= sample_years(t) && \
            lakes(lake_inds(i)).discovered == 0 && \
            runif(0,1) < pr_sample)
         {
            lakes(lake_inds(i)).discovered = sample_years(t);
         }
         if(t_vec(lake_inds(i)) > sample_years(t) && runif(0,1) < pr_sample)
         {
            lakes(lake_inds(i)).last_abs = sample_years(t);
         }
      }
   }
   
   //for(int i=1;i<=n_lakes;i++)
   //   cout << t_vec(i) << "\t" << lakes(i).last_abs << "\t" << lakes(i).discovered << "\n";
}

// Given the parameter names defined in the init_sims.ini file,
// produce a vector the current global values of those
// parameters (those used to generate the sim).
_vbc_vec<float> link_gen_pars(_vbc_vec<string> parname_str)
{
    _vbc_vec<float> parval(1,parname_str.UBound());
    for(int i = 1; i<= parname_str.UBound(); i++)
    {
        if(parname_str(i) == "alpha")
            parval(i) = log(glb_alpha);
        else if(parname_str(i) == "d")
        {
            parval(i) = d_par;
            traf_mat_precompute = TRUE;
            glb_traf_array_ind = lookup_traf_mat_array_index();
        }
        else if(parname_str(i) == "e")
        {
            parval(i) = e_par;
            traf_mat_precompute = TRUE;
            glb_traf_array_ind = lookup_traf_mat_array_index();
        }
        else if(parname_str(i) == "c")
            parval(i) = c_par;
        else if(parname_str(i) == "gamma")
            parval(i) = gamma_par;
        else if(parname_str(i) == "pdet")
        {
            parval(i) = pdet;
            fit_pdet = TRUE;
        }
        else if(parname_str(i) == "pdet_var0")
        {
            est_pdet = 0;
            parval(i) = pdet_pars(1);
            fit_pdet = TRUE;
            hetero_pdet = TRUE;
        }
        else if(parname_str(i) == "pdet_var1")
        {
            est_pdet = 1;
            parval(i) = pdet_pars(2);
            fit_pdet = TRUE;
            hetero_pdet = TRUE;
        }
        else if(parname_str(i) == "env0") //intercept
        {
            parval(i) = chem_pars(1);
            env = TRUE;
        }
        else if(parname_str(i) == "env1")
        {
            est_env = 1;
            parval(i) = chem_pars(2);
        }
        else if(parname_str(i) == "env2")
        {
            est_env = 2;
            parval(i) = chem_pars(3);
        }
        else if(parname_str(i) == "env3")
        {
            est_env = 3;
            parval(i) = chem_pars(4);
        }
        else if(parname_str(i) == "env4")
        {
            est_env = 4;
            parval(i) = chem_pars(5);
        }
        else if(parname_str(i) == "env5")
        {
            est_env = 5;
            parval(i) = chem_pars(6);
        }
        else if(parname_str(i) == "env6")
        {
            est_env = 6;
            parval(i) = chem_pars(7);
        }
        else if(parname_str(i) == "env7")
        {
            est_env = 7;
            parval(i) = chem_pars(8);
        }
        else if(parname_str(i) == "env8")
        {
            est_env = 8;
            parval(i) = chem_pars(9);
        }
        else if(parname_str(i) == "env9")
        {
            est_env = 9;
            parval(i) = chem_pars(10);
        }
        else if(parname_str(i) == "env10")
        {
            est_env = 10;
            parval(i) = chem_pars(11);
        }
        else if(parname_str(i) == "env11")
        {
            est_env = 11;
            parval(i) = chem_pars(12);
        }
        else if(parname_str(i) == "env12")
        {
            est_env = 12;
            parval(i) = chem_pars(13);
        }
        else if(parname_str(i) == "env13")
        {
            est_env = 13;
            parval(i) = chem_pars(14);
        }
        else if(parname_str(i) == "env14")
        {
            est_env = 14;
            parval(i) = chem_pars(15);
        }
        else if(parname_str(i) == "env15")
        {
            est_env = 15;
            parval(i) = chem_pars(16);
        }
        else if(parname_str(i) == "env16")
        {
            est_env = 16;
            parval(i) = chem_pars(17);
        }
        else if(parname_str(i) == "env17")
        {
            est_env = 17;
            parval(i) = chem_pars(18);
        }
        else if(parname_str(i) == "env18")
        {
            est_env = 18;
            parval(i) = chem_pars(19);
        }
        else if(parname_str(i) == "env19")
        {
            est_env = 19;
            parval(i) = chem_pars(20);
        }
        else if(parname_str(i) == "env20")
        {
            est_env = 20;
            parval(i) = chem_pars(21);
        }
        else if(parname_str(i) == "env1_2"){
            alpha_order = 2;
            parval(i) = chem_pars(est_env+2);
        }
        else if(parname_str(i) == "env2_2"){
            alpha_order = 2;
            parval(i) = chem_pars(est_env+3);
        }
        else if(parname_str(i) == "env3_2"){
            alpha_order = 2;
            parval(i) = chem_pars(est_env+4);
        }
        else if(parname_str(i) == "env4_2"){
            alpha_order = 2;
            parval(i) = chem_pars(est_env+5);
        }
        else if(parname_str(i) == "env5_2"){
            alpha_order = 2;
            parval(i) = chem_pars(est_env+6);
        }
        else if(parname_str(i) == "env6_2"){
            alpha_order = 2;
            parval(i) = chem_pars(est_env+7);
        }
        else if(parname_str(i) == "env7_2"){
            alpha_order = 2;
            parval(i) = chem_pars(est_env+8);
        }
        else if(parname_str(i) == "env8_2"){
            alpha_order = 2;
            parval(i) = chem_pars(est_env+9);
        }

        else if(parname_str(i) == "dummy")
            parval(i) = 0;
        else
        {
            cout << "Malformed parameter name\n"; 
            exit(1);
        }
    }
    return parval;

}


// Returns an n_lakes by time matrix of the probability of invasion for the current paramters and sim.
_vbc_vec<float> prob_inv()
{   
    _vbc_vec<float> pi(1,n_lakes,from_year+1,to_year);
    for(int i=1;i<=n_lakes;i++)
    {
        float alpha = calc_alpha(i);
        float log_not_inv = 0;
        for(int t=from_year+1;t<=to_year;t++)
        {
            update_pp_l_hood(i,t);
            log_not_inv += -pow(alpha*lakes(i).pp(t)+gamma_par,c_par);
            pi(i,t) = 1 - exp(log_not_inv);
        }
    }
    return pi;
}
