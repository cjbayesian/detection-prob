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
    // Sim lakes and sources //
    // <<< Instead just use 2EB
    inits();
    sim_inits();
    system("mkdir test_folder");
    // Read in Data // 
    cout << "Reading Spatial Data\n";
    read_data();

    // Set initial Params //
    cout << "Initiallizing...\n";
    // Initiallize (seed) the largest lake //
    // <<< Instead just use actual seed lakes
    // <<< But we will need to clear the current records of
    //     last_abs, last_obs, last_obs_uninv, discovered

    //for faster sims:
    n_lakes = 500;
    n_sources = 100;
    clear_traf_mat();
    int n_sim = 10000;

    one_time_traf = TRUE;

    for(int r=1;r<=100;r++)
    {
        string mf;
        string sf;
        stringstream out;
        out << r; 
        string s = out.str();
        mf = "test_folder/mean_field_" + s;
        sf = "test_folder/sim_" + s;
        ofstream mean_field_file(mf.c_str());
        ofstream sim_file(sf.c_str());

        for(int i=1;i<=n_lakes;i++)
        {
            //clear all but the seed lake(s)
            //lakes(i).Wj = log(lakes(i).Wj); //more normal size distr.
            if(lakes(i).discovered != from_year)
            {
                lakes(i).discovered = 0;
                lakes(i).invaded = 0;
                lakes(i).last_abs = 0;
                lakes(i).status_2010 = 0;
            }
        }
        init_state();
        init_t(); 
        calc_state();
        generate_params();

        // Simulate the invasion //
        calc_traf();
        calc_Uj();
        clear_t();
        calc_traf_mat();
        cout << "Simtest " << r << "\n";
        for(int s=1;s<=n_sim;s++)
        {
            sim_spread();
            for(int i=1;i<=n_lakes-1;i++)
                sim_file << t_vec(i) << ",";
            sim_file << t_vec(n_lakes) << "\n";
        }             
        cout << "Mean Field " << r << "\n";
        mean_field();
        cout << "Mean Field Done " << r << "\n";
        for(int t=from_year+1;t<=to_year;t++)
        {   
            for(int i=1;i<=n_lakes-1;i++)     
                mean_field_file << lakes(i).PInv(t) << ",";
            mean_field_file << lakes(n_lakes).PInv(t) << "\n";
        }
        
        mean_field_file.close();
        sim_file.close();
    }

    //////////////////////////////////////////////////

    return 0;
}


