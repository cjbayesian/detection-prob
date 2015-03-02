#include"grav.h"
#include"sim_grav.h"

int main(int argc,char *argv[])
{
    //string filename("b.mcmc");
    //cout << filename << " has " << wc_column(filename) << " columns.\n"; 

    args(argc,argv);
    if(output_folder != "")
        int status = mkdir(output_folder.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);


    // Sorting Example //
    _vbc_vec<float> nums(1,10);
    for(int i=1;i<=10;i++)
        nums(i) = runif(1,10);

    _vbc_vec<float> s_nums(1,10);
    _vbc_vec<int> rankO(1,10);
    _vbc_vec<int> rankN(1,10);
    sort(&nums,&s_nums,&rankO,&rankN);

    for(int i=1;i<=10;i++)
    {
        cout << fixed<< nums(i) << "\t" << s_nums(i) << "\t" << rankO(i) << "\t" << rankN(i) << "\t" << nums(rankO(i)) << "\n";
    }

    // -------------  //


    //string deep_folder("this/is/the/path/to/enlightenment");
    //int rrr = mkdir(deep_folder.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    //cout << "Deep Folder Created: " << rrr << "\n";
    // DATA HANDLING //
    inits();
    to_year = 2023;
    sim_inits();
    read_po_data();
    normalize_vars();
    n_sources = 100;
    clear_traf_mat();
    /*
    for(int i=1;i<=n_lakes;i++)
    {
        for(int j = 1;j<=n_chem_var-1;j++)
            cout << lakes(i).chem(j) << ",";
        cout <<lakes(i).chem(n_chem_var) << "\n";
    }
    */
    _vbc_vec<float> params1;
    parse_par_seeds(parseeds_file,&params1);
    link_pars(parnames,params1);
    for(int i=1;i<=n_lakes;i++)
    {
        cout << calc_alpha(i) << "\n";
    }
    //calc_traf_mat_array();
    //write_traf_mat_array();
    /*
    ofstream tmats("tmats_big.csv");
    for(int n =1;n<=tmat_res*tmat_res;n++)
    {
    for(int i=1;i<=n_lakes;i++)
    {
        for(int j=1;j<=n_lakes;j++) //start at 0 since the diagonal of traf_mat is not needed.
        {
            tmats << traf_mat_array(n,i,j) << ",";
        }
        tmats << n << "\n";
    }
        cout << n << "\n";
    }
    */
//    cerr << lakes(i).id << ": " << lakes(i).discovered << "\n";
//    cerr << "\n";
    // --------- //

    return 0;
}


