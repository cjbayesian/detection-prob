#ifndef MAN_GRAV_H
#define MAN_GRAV_H

////////////////////////////////////////////////////////////////////////////////
/////////////////////  Management Scenarios  ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Globals //


#endif;


// Set the management vectors //
_vbc_vec<int> set_top_n_source_lakes(int n, float theta, float phi, _vbc_vec<float> p_inv)
{
    //outflows is the net traffic from site i to
    // all other sites times probability that site i is invaded//
    _vbc_vec<float> outflows(1,n_lakes);
    for(int i=1;i<=n_lakes;i++)
    {
        // make sure no old policies are still in place //
        theta_vec(i) = 1;
        phi_vec(i) = 1;
        // ---------------------------------------------//

        outflows(i) = 0;
        for(int j=1;j<=n_lakes;j++)
        {   
            if(i != j)
            {
                for(int r=1;r<=tmat_res*tmat_res;r++)
                    outflows(i) += traf_mat_array(1,r,i,j);
            }
        }
    }

    // Sort the outflows (increasing order) //
    _vbc_vec<float> s_outflows(1,n_lakes);
    _vbc_vec<int> rankO(1,n_lakes);
    _vbc_vec<int> rankN(1,n_lakes);
    sort(&outflows,&s_outflows,&rankO,&rankN);


    // Impliment management at top n //    
    _vbc_vec<int> topn_index(1,n);
    for(int i=n_lakes;i >= n_lakes-(n-1); i--)
    {
        theta_vec(rankO(i)) = theta;
        phi_vec(rankO(i)) = phi;
        topn_index((n_lakes+1) - i) = rankO(i);
    }    
    // ---------------------------- //

    return topn_index;
}
