* Encoding: UTF-8.

****************************************************************************************
    *** CLUSTER-CIRC SIMU: cc_simu
****************************************************************************************
    *** Simulates data with perfect circumplex clusters in the population
    *** using parameters from the empirical data (sample size, number of clusters,
    *** number of items, empirical within-cluster range, mean item communality).
    *** Performs Cluster-Circ on the simulated data for comparison with results
    *** from Cluster-Circ Data. cc_simu can only be used after performing cc_data.


********************************************************************************************************************
    *** INSERTIONS: FILL IN SPECIFICATIONS OF YOUR DATA: DEMONSTRATION ON EXEMPLARY DATA
********************************************************************************************************************
    *** Replace each element with "$INSERT" with the relevant information from your data
    
    *** $INSERT_path:            Working directory/path with the dataset and where ClusterCirc results should be saved. Same as in cc_data (line 26).
    *** $INSERT_n:                  Number of subjects in the sample (line 42)
    *** $INSERT_samples:    Number of samples for the simulation. Recommendation: 100-500 (depending on computing capacity).
    ***                                          Must be inserted two times (lines 36, 43).
    *** $INSERT_m+2:            Number of items in the data + 2 (e.g. 20 if your data has 18 items). Needs to be inserted 3 times (lines 46; 63; 70)
    
    *** Here: Only $INSERT_path needs to be replaced. Otherwise, parameters for the exemplary data (n = 300,  18 variables -> m+2 = 20 ) and samples = 500 are used.

cd "$INSERT_path".

*****************************************
    *** CREATE POPULATION DATA
*****************************************
 
    * Normally distributed z-values for m error factors and 2 circumplex factors.
    * Perform PCA to orthogonalize factors and use regression factor scores (U,F).
 
get file = "help_CC_data.sav".
compute samples = 500.
save outfile = "help_CC_data.sav".

SET RNG=MT MTINDEX=3391001000 workspace=200000 MXLOOPS=300000000.

INPUT PROGRAM.
    compute n = 300.
    compute samples = 500.
    compute n_pop = n*samples.
    LOOP #1 = 1 to n_pop.
        DO REPEAT z = z1 to z20.
            COMPUTE z = NORMAL(1).
        END REPEAT.
        END CASE.
    END LOOP.
    END FILE.
END INPUT PROGRAM.
EXECUTE.

save outfile = "z_normal.sav"
    /drop  n samples n_pop.

get file = "z_normal.sav".

OMS /SELECT ALL /DESTINATION VIEWER=NO.
FACTOR /VARIABLES ALL
    /PRINT INITIAL EXTRACTION
    /CRITERIA FACTORS(20) ITERATE(500)
    /EXTRACTION PC
    /ROTATION NOROTATE
    /SAVE REG(ALL x)
    /METHOD=CORRELATION.
OMSEND.

delete variables z1 to z20.
save outfile ="Pop_orthfactors.sav".
get file = "Pop_orthfactors.sav".


**************************************
    **** END OF INSERTIONS
**************************************

set mxloops = 300000000.

MATRIX.

    get factors /file= "Pop_orthfactors.sav"
        /missing omit.
    get p /variables p /file = "help.sav".
    get m /variables m /file = "help.sav".
    get n /variables n /file = "help.sav".
    get w_range /variables wrange_c /file =  "help_CC_data.sav".
    get h_sq /variables Hsq_mn /file = "help_CC_data.sav".
    get n_simu /variables samples /file = "help_CC_data.sav".
    
    compute c_wrange = make(p, 1, w_range).
    compute h = sqrt(make(m, 1, h_sq)).
    compute k = 2.
    compute n_pop = n*n_simu.

    * The number of items per cluster mc in the population is chosen to be roughly 
    * equal in all clusters. If m/p is not a natural number, some clusters have fewer 
    * and some have more items:
    * mp_s = small mp, mp_l = large mp,  dec = decimals (for large and small mp)
    * mc = actual number of items in cluster c
    * p1 = number of clusters with small mc, p2 = number of clusters with large mc

    compute mp = m/p.
    compute mp_s = trunc(mp).
    compute mp_l = mp_s + 1.
    compute dec_l = mp - mp_s.

    do if dec_l = 0.
        compute mc = make(p,1,mp).
    else if dec_l > 0.
        compute dec_s = 1 - dec_l.
        compute p1 = p* dec_s.
        compute p2 = p*dec_l.
        compute mc_s = make(p1,1,mp_s).
        compute mc_l = make(p2,1,mp_l).
        compute mc = {mc_s;mc_l}.
    end if.

    compute U = make(n_pop, m, 0).
    compute F = make(n_pop,2,0).
    loop j = 1 to m.
        compute U(:,j) = factors(:,j).
    end loop.
    loop j = 1 to 2.
        compute F(:,j) = factors(:,m+j).
    end loop.

     * Introduce circumplexity: Translate perfect circumplex angles into population loadings.
     * Use reg factor scores and make data with var = F*T(A) + U*T(D).

    compute pi = 4*ARTAN(1).
    compute cl_dis = 360/p.
    compute ones = make (p,1,1).
    compute w_dis = c_wrange&/(mc-ones).

    loop c = 1 to p.
        compute th = make(mc(c),1,0).
        loop i = 1 to mc(c).
            compute th(i) = cl_dis*(c-1) + (i-(mc(c)+1)/2)*w_dis(c).
        end loop.
        do if c = 1.
            compute theta_d = th.
        else if c > 1.
            compute theta_d = {theta_d; th}.
        end if.
    end loop.

    loop i = 1 to m.
        do if theta_d(i) < 0.
            compute theta_d(i) = 360 + theta_d(i).
        end if.
    end loop.

    compute theta = theta_d*(pi/180).

    * Angles > 2*pi (>360°) need to be converted

    loop i = 1 to m.
        do if theta(i) > 2*pi.
            compute theta(i) = theta(i) - 2*pi.
        end if.
    end loop.

    * Conversion with sin/cos works only for theta >= 0 and <= pi/2 (90°).
    * Help vector theta_h for computation of population loadings.

    compute theta_h = theta.
    loop i = 1 to m.
        do if theta(i) >pi/2 and theta(i) <= pi.
            compute theta_h(i) = pi - theta(i).
        end if.
        do if theta(i) >pi and theta(i) <= 3/2*pi.
            compute theta_h(i) = theta(i) - pi.
        end if.
        do if theta(i) >3/2*pi and theta(i) <= 2*pi.
            compute theta_h(i) = 2*pi-theta(i).
        end if.
    end loop.

    compute A = make(m,2,0).
    loop i = 1 to m.
        do if theta(i) >= 0 and theta(i) <=pi/2.
            compute A(i,1) = h(i)*sin(theta_h(i)).
            compute A(i,2) = h(i)*cos(theta_h(i)).
        end if.
        do if theta(i) > pi/2 and theta(i) <=pi.
            compute A(i,1) = h(i)*sin(theta_h(i)).
            compute A(i,2) = -(h(i)*cos(theta_h(i))).
        end if.
        do if theta(i) > pi and theta(i) <=3/2*pi.
            compute A(i,1) = -(h(i)*sin(theta_h(i))).
            compute A(i,2) = -(h(i)*cos(theta_h(i))).
        end if.
        do if theta(i) >3/2*pi and theta(i) <=2*pi.
            compute A(i,1) = -(h(i)*sin(theta_h(i))).
            compute A(i,2) = h(i)*cos(theta_h(i)).
        end if.
    end loop.

    * Get error loadings D for var = F*T(A) + U*T(D)

    compute ones = make(1,m,1).
    compute ones = mdiag(ones).
    compute AtA = A*T(A).
    compute AtA = mdiag(diag(AtA)).
    compute D = sqrt(ones-AtA).

    * Compute variables (simulated values for each subject in the population):

    compute Var = F*T(A) + U*T(D).

    loop i = 1 to n_simu.
        do if i = 1.
            compute sample= make(n, 1, i).
        else if i > 1.
            compute sample = {sample; make(n, 1, i)}.
        end if.
    end loop.

    save mc /outfile = "Pop_mc.sav".
    save Var /outfile = "simudata.sav".
    save {sample, Var} /outfile ="simudata_sample.sav"
        /variables sample var00001 to var10000.

END MATRIX.


****************************************
    *** PCA ON POPULATION DATA
****************************************

get file = "simudata.sav".

OMS /SELECT ALL /DESTINATION VIEWER=NO.
FACTOR
    /VARIABLES ALL
    /MISSING LISTWISE
    /MATRIX OUT(fac="Pop_A.sav")
    /ANALYSIS ALL
    /PRINT INITIAL EXTRACTION ROTATION
    /CRITERIA FACTORS(2) ITERATE(250)
    /EXTRACTION PC
    /CRITERIA ITERATE(500)
    /ROTATION NOROTATE
    /METHOD=CORRELATION.
OMSEND.

get file = "Pop_A.sav".
delete variables ROWTYPE_ FACTOR_.
flip variables all.
rename variables var001 = F1 var002 = F2.
save outfile = "Pop_A.sav".
get file =  "Pop_A.sav".

***************************************
    *** PCA ON SAMPLES
***************************************

OMS /SELECT ALL /DESTINATION VIEWER=NO.
get file= "simudata_sample.sav".
spssinc select variables MACRONAME = "!allvar"
    /PROPERTIES PATTERN="var*".

SORT CASES  BY sample.
SPLIT FILE SEPARATE BY sample.
FACTOR  /VARIABLES !allvar
    /MATRIX=OUT(FAC="Simu_A.sav")
    /MISSING LISTWISE   /PRINT KMO  /CRITERIA FACTORS(2) ITERATE(250)
    /EXTRACTION PC  /CRITERIA ITERATE(250) NOKAISER /ROTATION NOROTATE /METHOD CORRELATION.

OMSEND.
SPLIT FILE OFF.

get file = "Simu_A.sav".
delete variables ROWTYPE_ FACTOR_ sample.
flip variables all.
delete variables CASE_LBL.
save outfile = "Simu_A.sav".
get file "Simu_A.sav".

output close *.


*******************************************************************
    *** CLUSTER-CIRC ON POPULATION LOADINGS
*******************************************************************

set mxloops = 300000000.
set mdisplay TABLES.

MATRIX.

    get A /variables F1 F2 /file = "Pop_A.sav".
    get p /variables p /file = "help.sav".
    get m /variables m /file = "help.sav".
    get q /variables q /file = "help.sav".
    get n /variables n /file = "help.sav".
    get n_simu /variables samples /file = "help_CC_data.sav".
    get mc_id /file = "Pop_mc.sav".

  *************************
        *** PREPARATION
    *************************
        * Kaiser-normalization of loadings for simpler computation of angles

    compute h_sq = mdiag(rssq(A)).
    compute h_rt = sqrt(h_sq).
    compute h_sqv = rsum(h_sq).
    compute h_rtv = rsum(h_rt).
    compute hsq_mn = msum(h_sq)*(1/m).
    compute A_k = inv(h_rt)*A.

        * Compute and sort theta: item angles in degrees
        * r = radians, p = positive
        
    compute pi=4*ARTAN(1).
    compute A_pos = abs(A_k).
    compute th_rp = arsin(A_pos(1:m,1)).
    compute th_r = make(m,1,0).
    compute theta = make(m,1,0).

        * Computation of angles depends on the quadrant (loadings positive/negative)

    loop i = 1 to m.
        do if A_k(i,1) >= 0 and A_k(i,2) >= 0.
            compute th_r (i) = th_rp(i).
        end if.
        do if A_k(i,1) >= 0 and A_k(i,2) < 0.
            compute th_r(i) = pi - th_rp(i).
        end if.
        do if A_k(i,1) < 0 and A_k(i,2) < 0.
            compute th_r(i) = pi + th_rp(i).
        end if.
        do if A_k(i,1) < 0 and A_k(i,2) >= 0.
            compute th_r(i) = 2*pi - th_rp(i).
        end if.
    compute theta(i) = th_r(i)*180/pi. 
    end loop.

    * Sort theta and keep item number in ival (for later re-assignment)

    compute rk_th = grade(theta).
    compute ival = make(m,3,0).

    loop i1 = 1 to m.
        loop i2 = 1 to m.
            do if rk_th(i1) = i2.
                compute ival(i2,2) = theta(i1).
                compute ival(i2,1) = i1.
                compute ival(i2,3) = h_sqv(i1).
            end if.
        end loop.
    end loop.

    ***********************************************
        *** CLUSTER-CIRC ALGORITHM
    ***********************************************
    
    compute spacingh = 361.

     LOOP d = 0 to rnd(360*q/p).

        compute ci_h = make(m,1,0).
        compute c_m = make(p,1,0).
        compute c_no = make(p,1,0).
        compute c_min = make(p,1,0).
        compute c_max = make(p,1,0).
        compute c_rng = make(p,1,0).
        compute c_ang = make(p,1,0).

        * Check if item falls within the range of a cluster. Adjust for max. 360°.

        loop c = 1 to p.
            loop i = 1 to m.

                compute thmin = (c-1)*360/p +d/q.
                compute thmax = c*360/p + d/q.

                do if thmin <= 360 and thmax<= 360.
                    do if ival(i,2) >= thmin and ival(i,2)< thmax.
                        compute ci_h(i) = c.
                        compute c_m(c) = c_m(c) + 1.
                    end if.
                end if.

                do if thmin <=360 and thmax > 360.
                    compute thmax_n = thmax-360.
                    do if ival(i,2) >=thmin and ival(i,2) <=360 or ival(i,2) >=0 and ival(i,2)<thmax_n.
                        compute ci_h(i) = c.
                        compute c_m(c) = c_m(c) + 1.
                    end if.
                end if.

                do if thmin > 360 and thmax > 360.
                    compute thmin_n = thmin-360.
                    compute thmax_n = thmax-360.
                    do if ival(i,2) >= thmin_n and ival(i,2)< thmax_n.
                        compute ci_h(i) = c.
                        compute c_m(c) = c_m(c) + 1.
                    end if.
                end if.

            end loop.
        end loop.

        * If c_m = 0: One cluster is empty. Prepare dismissal of solution.

        loop c = 1 to p.
            do if c_m(c) = 0.
                compute c_m(c) = 99.
            end if.
        end loop.

        * Compute cluster angle as the center between the outer items in cluster

        loop c = 1 to p.

            compute c_min(c) = 361.
            compute c_max(c) = 0.

            loop i = 1 to m.
                do if ci_h(i) = c and ival(i,2) <= c_min(c).
                    compute c_min(c) = ival(i,2).
                end if.
                do if ci_h(i) = c and ival(i,2) >= c_max(c).
                    compute c_max(c) = ival(i,2).
                end if.
            end loop.

            do if c_max(c)-c_min(c) < 180.
                compute c_ang(c) = (c_max(c)+c_min(c))/2.
                compute c_rng(c) = c_max(c)-c_min(c).
            end if.
      
            * Special case: Cluster at approx. 0° could have a range of > 180.
            * Change item angles with help objects.

            do if c_max(c)-c_min(c) > 180.

                compute ival_h = ival(:,2).
                compute c_minh = 361.
                compute c_maxh = 0.

                loop i = 1 to m.
                    do if ival(i,2) > 180.
                        compute ival_h(i) = ival(i,2)-360.
                    end if.
                    do if ci_h(i) = c and ival_h(i) <= c_minh.
                        compute c_minh = ival_h(i).
                    end if.
                    do if ci_h(i) = c and ival_h(i) >= c_maxh.
                        compute c_maxh = ival_h(i).
                    end if.
                end loop.

                compute c_max(c) = c_maxh.
                compute c_min(c) = 360+c_minh.
                compute c_rng(c) = c_maxh -(c_minh).
                compute c_ang(c) = (c_minh + c_maxh)/2.

            end if.

        end loop.
    
        * Clusters and items need to be sorted before computing spacing indices

        compute c_rnk = grade(c_ang).
        compute cval = make(p,4,0).

        loop c1 = 1 to p.
            loop c2 = 1 to p.
                do if c_rnk(c2) = c1.
                    compute cval(c1,1) = c2.
                    compute cval(c1,2) = c_m(c2).
                    compute cval(c1,3) = c_ang(c2).
                    compute cval(c1,4) = c_rng(c2).
                end if.
            end loop.
        end loop.

        compute c_i = make(m,1,0).
        compute cvalh = cval.

        loop i = 1 to m.
            loop c  = 1 to p.
                do if ci_h(i) = cval(c,1).
                    compute c_i(i) = c.
                    compute cvalh(c,1) = c.
                end if.
            end loop.
        end loop.

        * Help vectors for computation of distances:
        * Allow for negative angles (n) in first cluster

        compute ival_n = {c_i, ival}.
        loop i = 1 to m.
            do if c_i(i) = 1 and ival(i,2) > 180.
                compute ival_n(i,3) = ival_n(i,3) - 360.
            end if.
        end loop.

        compute rk_iang = grade(ival_n(:,3)).
        compute ival_h = make(m,4,0).

        loop i1 = 1 to m.
            loop i2 = 1 to m.
                do if rk_iang(i2) = i1.
                    compute ival_h(i1,:) = ival_n(i2,:).
                end if.
            end loop.
        end loop.

        * Compute spacing_h for each division

        compute ic_dis = make(m,p,0).
        compute space = 360/p.
        compute ic_dev = make(m,p,0).
        compute ic_devp = make(m,p,0).
        compute ic_dcom = make(m,p,0).

        loop i = 1 to m.
            loop c1 = 1 to p.
                compute c2 = ival_h(i,1).
                compute i_ang =  ival_h(i,3).
                compute c_ang = cvalh(c1,3).
                compute i_com = ival_h(i,4).
                compute ic_dis(i,c1) = i_ang - c_ang.
                compute id_dis = (c2-1)*space - (c1-1)*space.
                compute ic_dev(i,c1) = ic_dis(i,c1) - id_dis.
                compute ic_devp(i,c1) = ic_dev(i,c1)/space.
                compute ic_dcom(i,c1) = ic_devp(i,c1)* sqrt(i_com).
            end loop.
        end loop.

        * Item-spacing:

        compute ispc_sq = rssq(ic_devp)/p.
        compute ispc = sqrt(ispc_sq).

        * With communalities (h) for spacing index (not interpretable on item level)

        compute ispc_hsq = rssq(ic_dcom)/(p*hsq_mn).
        compute ispc_h = sqrt(ispc_hsq).

        * Overall spacing

        compute spc_sq  = csum(ispc_sq)/m.
        compute spc = sqrt(spc_sq).

        compute spc_hsq = csum(ispc_hsq)/m.
        compute spc_h = sqrt(spc_hsq).

        * Dismiss partitions with empty clusters by making spc_h larger than
        * the initial spacingh (361). If this happens for all possible divisions,
        * the number of clusters is too large.

        loop c = 1 to p.
            do if c_m(c) = 99.
                compute spc_h = 50000.
            end if.
        end loop.

        do if spc_h < spacingh.
            compute spacingh = spc_h.
            compute spacing  = spc.
            compute items = {ival_h, ispc}.
            compute clusters = cvalh.
            compute ic_dist = ic_dis.
        end if.

    END LOOP.

    do if spacingh = 361.
        print/Title "Cluster-Circ could not finish, at least one of the clusters is empty. Try a smaller number of clusters or include more variables.".

    else if spacingh < 361.

    * Between-cluster spacing

        compute c_dis = make(p,p,0).
        compute space = 360/p.
        compute c_dev = make(p,p,0).
        compute c_devp = make(p,p,0).

        loop c1 = 1 to p.
            loop c2 = 1 to p.
                compute c_dis(c1,c2) = clusters(c1,3)-clusters(c2,3).
                compute id_dis = (c1-1)*space - (c2-1)*space.
                compute c_dev(c1,c2) = (c_dis(c1,c2) - id_dis).
                compute c_devp(c1,c2) = c_dev(c1,c2)/space.
            end loop.
        end loop.

        compute cbcs_psq = rssq(c_devp)/(p-1).
        compute cbcs_p = sqrt(cbcs_psq).

        compute clusters = {clusters, cbcs_p}.

        * Overall between-cluster spacing

        compute bcs_p = sqrt(csum(cbcs_psq)/p).

        * Within-cluster proximity

        compute ic_disw = make(m,1,0).
        compute ic_disp = make(m,1,0).

        compute cwcp_sqp = make(p,1,0).
        compute cwcp_p = make(p,1,0).

        loop c = 1 to p.

            compute c_m = clusters(c,2).
            compute c_icdisp = make(m,1,0).

            loop i = 1 to m.
                do if items(i,1) = c.
                    compute ic_disw(i) = abs(ic_dist(i,c)).
                    compute ic_disp(i) = (ic_disw(i)/(360/(2*p))).
                    compute c_icdisp(i) = ic_disp(i).
                end if.
            end loop.

            compute cwcp_sqp(c) = cssq(c_icdisp)/c_m.
            compute cwcp_p(c) = sqrt(cwcp_sqp(c)).

        end loop.

        * Overall within cluster spacing

        compute wcp_p = sqrt(csum(cwcp_sqp)/p).

        * Final values

        compute clusters= {clusters, cwcp_p}.
        compute items ={items, ic_disp}.

        * Recompute negative angles

        loop i = 1 to m.
            do if items(i,3) < 0.
                compute items(i,3) = 360+items(i,3).
            end if.
        end loop.

        loop c = 1 to p.
            do if clusters(c,3) < 0.
                compute clusters(c,3) = 360+clusters(c,3).
            end if.
        end loop.

        * Final overall Cluster-Circ indices
        
        compute overall = {spacingh, spacing,  bcs_p, wcp_p}.

        ******************************
            *** END CLUSTER-CIRC
        ******************************

        ********************************
            *** SORTING CORRECT?
        ********************************
            * clus_id = ideal cluster sorting

        loop c = 1 to p.
            compute c_id = make(mc_id(c),1,c).
            do if c = 1.
                compute clus_id = c_id.
            else if c > 1.
                compute clus_id = {clus_id;c_id}.
            end if.
        end loop.

            * Each entry in 'sort_it' indicates whether items are clustered together correctly (1) or wrongly (0).
            * A single zero entry indicates that the intended structure has not been found.
            * Position in clus_id, items(i,2)  = item number.
            * items(i,1) is the cluster to which the item has been sorted.

        compute sort_it = make(m,m,0).

        loop i1 = 1 to m.
            loop i2 = 1 to m.
                do if i1 <> i2.
                    do if clus_id(items(i1,2)) = clus_id(items(i2,2)) and items(i1,1) = items(i2,1).
                        compute sort_it(i1,i2) = 1.
                    end if.
                    do if clus_id(items(i1,2)) = clus_id(items(i2,2)) and items(i1,1) <> items(i2,1).
                        compute sort_it(i1,i2) = 0.
                    end if.
                    do if clus_id(items(i1,2)) <> clus_id(items(i2,2)) and items(i1,1) = items(i2,1).
                        compute sort_it(i1,i2) = 0.
                    end if.
                    do if clus_id(items(i1,2)) <> clus_id(items(i2,2)) and items(i1,1) <> items(i2,1).
                        compute sort_it(i1,i2) = 1.
                    end if.
                end if.
            end loop.
        end loop.

        do if csum(rsum(sort_it)) = m*(m-1).
            compute sortcorr = 1.
        else if csum(rsum(sort_it)) <m*(m-1).
            compute sortcorr = 0.
        end if.

        compute overall = {overall}.

        compute simu = {n_simu, n, m, p, q}.

        print /Title = "RESULTS CLUSTER-CIRC SIMU".
        print simu 
            /rlabels = "Parameters of the simulation"
           /clabels = "Number of simulated samples", "Sample size", "Items", "Clusters", "Precision index".

        print /Title = "POPULATION (SIMULATED)".

        print /Title = "The following tables show results for a simulated population with perfect circumplex spacing of clusters,".
        print /Title = "adapted for the specifications of the empirical data: Sample size, number of clusters, number of items,".
        print /Title = "empirical within-cluster range, mean item communality).".

        print overall
            /rlabels = " Population"
            /clabels = "Spacing (with h²)", "Spacing", "Between-cluster spacing",  "Within-cluster proximity ".

        do if sortcorr = 1.
            print /Title = "Cluster-Circ found the intended circumplex clusters in the population.".
        else if sortcorr = 0.
            print /Title = "Cluster-Circ did not find the intended circumplex clusters in the population".
        end if.
    
        print clusters
            /rlabels = " Population"
            /clabels = "Cluster", "Items", "Angles", "Range",  "Between-cluster spacing ",  "Within-cluster proximity".

        print items
            /rlabels = " Population"
            /clabels = "Cluster", "Item", "Angle", "Communality",  "Item-cluster spacing ", "Distance to cluster center (0-1)".

    end if.


END MATRIX.


**************************************
    ** SAMPLES: CLUSTER-CIRC
**************************************

get file = "Pop_A.sav".

set mxloops = 3000000000000.

MATRIX.

    get A_Pop /variables F1 F2 /file = "Pop_A.sav".
    get A_all /file "Simu_A.sav".
    get spch_dat /variables spacingh /file = "help_CC_data.sav".
    get n_simu /variables samples /file = "help_CC_data.sav".
    get p /variables p /file = "help.sav".
    get m /variables m /file = "help.sav".
    get n /variables n /file = "help.sav".
    get q /variables q /file = "help.sav".
    get mc_id /file = "Pop_mc.sav".

    compute k = 2.
    compute n_pop = n*n_simu.

    LOOP i_simu = 1 to n_simu.

        compute A = make(m,k,0).

        loop g = 1 to m.
            loop j = 1 to k.
                compute i_sample= (i_simu-1)*k+j.
                compute A(g, j)= A_all(g, i_sample).
            end loop.
        end loop.

 
  *************************
        *** PREPARATION
    *************************
        * Kaiser-normalization of loadings for simpler computation of angles

    compute h_sq = mdiag(rssq(A)).
    compute h_rt = sqrt(h_sq).
    compute h_sqv = rsum(h_sq).
    compute h_rtv = rsum(h_rt).
    compute hsq_mn = msum(h_sq)*(1/m).
    compute A_k = inv(h_rt)*A.

        * Compute and sort theta: item angles in degrees
        * r = radians, p = positive
        
    compute pi=4*ARTAN(1).
    compute A_pos = abs(A_k).
    compute th_rp = arsin(A_pos(1:m,1)).
    compute th_r = make(m,1,0).
    compute theta = make(m,1,0).

        * Computation of angles depends on the quadrant (loadings positive/negative)

    loop i = 1 to m.
        do if A_k(i,1) >= 0 and A_k(i,2) >= 0.
            compute th_r (i) = th_rp(i).
        end if.
        do if A_k(i,1) >= 0 and A_k(i,2) < 0.
            compute th_r(i) = pi - th_rp(i).
        end if.
        do if A_k(i,1) < 0 and A_k(i,2) < 0.
            compute th_r(i) = pi + th_rp(i).
        end if.
        do if A_k(i,1) < 0 and A_k(i,2) >= 0.
            compute th_r(i) = 2*pi - th_rp(i).
        end if.
    compute theta(i) = th_r(i)*180/pi. 
    end loop.

    * Sort theta and keep item number in ival (for later re-assignment)

    compute rk_th = grade(theta).
    compute ival = make(m,3,0).

    loop i1 = 1 to m.
        loop i2 = 1 to m.
            do if rk_th(i1) = i2.
                compute ival(i2,2) = theta(i1).
                compute ival(i2,1) = i1.
                compute ival(i2,3) = h_sqv(i1).
            end if.
        end loop.
    end loop.

    ***********************************************
        *** CLUSTER-CIRC ALGORITHM
    ***********************************************
    
    compute spacingh = 361.

     LOOP d = 0 to rnd(360*q/p).

        compute ci_h = make(m,1,0).
        compute c_m = make(p,1,0).
        compute c_no = make(p,1,0).
        compute c_min = make(p,1,0).
        compute c_max = make(p,1,0).
        compute c_rng = make(p,1,0).
        compute c_ang = make(p,1,0).

        * Check if item falls within the range of a cluster. Adjust for max. 360°.

        loop c = 1 to p.
            loop i = 1 to m.

                compute thmin = (c-1)*360/p +d/q.
                compute thmax = c*360/p + d/q.

                do if thmin <= 360 and thmax<= 360.
                    do if ival(i,2) >= thmin and ival(i,2)< thmax.
                        compute ci_h(i) = c.
                        compute c_m(c) = c_m(c) + 1.
                    end if.
                end if.

                do if thmin <=360 and thmax > 360.
                    compute thmax_n = thmax-360.
                    do if ival(i,2) >=thmin and ival(i,2) <=360 or ival(i,2) >=0 and ival(i,2)<thmax_n.
                        compute ci_h(i) = c.
                        compute c_m(c) = c_m(c) + 1.
                    end if.
                end if.

                do if thmin > 360 and thmax > 360.
                    compute thmin_n = thmin-360.
                    compute thmax_n = thmax-360.
                    do if ival(i,2) >= thmin_n and ival(i,2)< thmax_n.
                        compute ci_h(i) = c.
                        compute c_m(c) = c_m(c) + 1.
                    end if.
                end if.

            end loop.
        end loop.

        * If c_m = 0: One cluster is empty. Prepare dismissal of solution.

        loop c = 1 to p.
            do if c_m(c) = 0.
                compute c_m(c) = 99.
            end if.
        end loop.

        * Compute cluster angle as the center between the outer items in cluster

        loop c = 1 to p.

            compute c_min(c) = 361.
            compute c_max(c) = 0.

            loop i = 1 to m.
                do if ci_h(i) = c and ival(i,2) <= c_min(c).
                    compute c_min(c) = ival(i,2).
                end if.
                do if ci_h(i) = c and ival(i,2) >= c_max(c).
                    compute c_max(c) = ival(i,2).
                end if.
            end loop.

            do if c_max(c)-c_min(c) < 180.
                compute c_ang(c) = (c_max(c)+c_min(c))/2.
                compute c_rng(c) = c_max(c)-c_min(c).
            end if.
      
            * Special case: Cluster at approx. 0° could have a range of > 180.
            * Change item angles with help objects.

            do if c_max(c)-c_min(c) > 180.

                compute ival_h = ival(:,2).
                compute c_minh = 361.
                compute c_maxh = 0.

                loop i = 1 to m.
                    do if ival(i,2) > 180.
                        compute ival_h(i) = ival(i,2)-360.
                    end if.
                    do if ci_h(i) = c and ival_h(i) <= c_minh.
                        compute c_minh = ival_h(i).
                    end if.
                    do if ci_h(i) = c and ival_h(i) >= c_maxh.
                        compute c_maxh = ival_h(i).
                    end if.
                end loop.

                compute c_max(c) = c_maxh.
                compute c_min(c) = 360+c_minh.
                compute c_rng(c) = c_maxh -(c_minh).
                compute c_ang(c) = (c_minh + c_maxh)/2.

            end if.

        end loop.
    
        * Clusters and items need to be sorted before computing spacing indices

        compute c_rnk = grade(c_ang).
        compute cval = make(p,4,0).

        loop c1 = 1 to p.
            loop c2 = 1 to p.
                do if c_rnk(c2) = c1.
                    compute cval(c1,1) = c2.
                    compute cval(c1,2) = c_m(c2).
                    compute cval(c1,3) = c_ang(c2).
                    compute cval(c1,4) = c_rng(c2).
                end if.
            end loop.
        end loop.

        compute c_i = make(m,1,0).
        compute cvalh = cval.

        loop i = 1 to m.
            loop c  = 1 to p.
                do if ci_h(i) = cval(c,1).
                    compute c_i(i) = c.
                    compute cvalh(c,1) = c.
                end if.
            end loop.
        end loop.

        * Help vectors for computation of distances:
        * Allow for negative angles (n) in first cluster

        compute ival_n = {c_i, ival}.
        loop i = 1 to m.
            do if c_i(i) = 1 and ival(i,2) > 180.
                compute ival_n(i,3) = ival_n(i,3) - 360.
            end if.
        end loop.

        compute rk_iang = grade(ival_n(:,3)).
        compute ival_h = make(m,4,0).

        loop i1 = 1 to m.
            loop i2 = 1 to m.
                do if rk_iang(i2) = i1.
                    compute ival_h(i1,:) = ival_n(i2,:).
                end if.
            end loop.
        end loop.

        * Compute spacing_h for each division

        compute ic_dis = make(m,p,0).
        compute space = 360/p.
        compute ic_dev = make(m,p,0).
        compute ic_devp = make(m,p,0).
        compute ic_dcom = make(m,p,0).

        loop i = 1 to m.
            loop c1 = 1 to p.
                compute c2 = ival_h(i,1).
                compute i_ang =  ival_h(i,3).
                compute c_ang = cvalh(c1,3).
                compute i_com = ival_h(i,4).
                compute ic_dis(i,c1) = i_ang - c_ang.
                compute id_dis = (c2-1)*space - (c1-1)*space.
                compute ic_dev(i,c1) = ic_dis(i,c1) - id_dis.
                compute ic_devp(i,c1) = ic_dev(i,c1)/space.
                compute ic_dcom(i,c1) = ic_devp(i,c1)* sqrt(i_com).
            end loop.
        end loop.

        * Item-spacing:

        compute ispc_sq = rssq(ic_devp)/p.
        compute ispc = sqrt(ispc_sq).

        * With communalities (h) for spacing index (not interpretable on item level)

        compute ispc_hsq = rssq(ic_dcom)/(p*hsq_mn).
        compute ispc_h = sqrt(ispc_hsq).

        * Overall spacing

        compute spc_sq  = csum(ispc_sq)/m.
        compute spc = sqrt(spc_sq).

        compute spc_hsq = csum(ispc_hsq)/m.
        compute spc_h = sqrt(spc_hsq).

        * Dismiss partitions with empty clusters by making spc_h larger than
        * the initial spacingh (361). If this happens for all possible divisions,
        * the number of clusters is too large.

        loop c = 1 to p.
            do if c_m(c) = 99.
                compute spc_h = 50000.
            end if.
        end loop.

        do if spc_h < spacingh.
            compute spacingh = spc_h.
            compute spacing  = spc.
            compute items = {ival_h, ispc}.
            compute clusters = cvalh.
            compute ic_dist = ic_dis.
        end if.

    END LOOP.

    do if spacingh = 361.
        print/Title "Cluster-Circ could not finish, at least one of the clusters is empty. Try a smaller number of clusters or include more variables.".

    else if spacingh < 361.

    * Between-cluster spacing

        compute c_dis = make(p,p,0).
        compute space = 360/p.
        compute c_dev = make(p,p,0).
        compute c_devp = make(p,p,0).

        loop c1 = 1 to p.
            loop c2 = 1 to p.
                compute c_dis(c1,c2) = clusters(c1,3)-clusters(c2,3).
                compute id_dis = (c1-1)*space - (c2-1)*space.
                compute c_dev(c1,c2) = (c_dis(c1,c2) - id_dis).
                compute c_devp(c1,c2) = c_dev(c1,c2)/space.
            end loop.
        end loop.

        compute cbcs_psq = rssq(c_devp)/(p-1).
        compute cbcs_p = sqrt(cbcs_psq).

        compute clusters = {clusters, cbcs_p}.

        * Overall between-cluster spacing

        compute bcs_p = sqrt(csum(cbcs_psq)/p).

        * Within-cluster proximity

        compute ic_disw = make(m,1,0).
        compute ic_disp = make(m,1,0).

        compute cwcp_sqp = make(p,1,0).
        compute cwcp_p = make(p,1,0).

        loop c = 1 to p.

            compute c_m = clusters(c,2).
            compute c_icdisp = make(m,1,0).

            loop i = 1 to m.
                do if items(i,1) = c.
                    compute ic_disw(i) = abs(ic_dist(i,c)).
                    compute ic_disp(i) = (ic_disw(i)/(360/(2*p))).
                    compute c_icdisp(i) = ic_disp(i).
                end if.
            end loop.

            compute cwcp_sqp(c) = cssq(c_icdisp)/c_m.
            compute cwcp_p(c) = sqrt(cwcp_sqp(c)).

        end loop.

        * Overall within cluster spacing

        compute wcp_p = sqrt(csum(cwcp_sqp)/p).

        * Final values

        compute clusters= {clusters, cwcp_p}.
        compute items ={items, ic_disp}.

        * Recompute negative angles

        loop i = 1 to m.
            do if items(i,3) < 0.
                compute items(i,3) = 360+items(i,3).
            end if.
        end loop.

        loop c = 1 to p.
            do if clusters(c,3) < 0.
                compute clusters(c,3) = 360+clusters(c,3).
            end if.
        end loop.

        * Final overall Cluster-Circ indices
        
        compute overall = {spacingh, spacing,  bcs_p, wcp_p}.


       *******************************************************
          *** OVERALL VALUES FOR SIMULATION
       ******************************************************
                
       ***************************************
          *** SORTING CORRECT?
       ***************************************

            loop c = 1 to p.
                compute c_id = make(mc_id(c),1,c).
                do if c = 1.
                    compute clus_id = c_id.
                else if c > 1.

                    compute clus_id = {clus_id;c_id}.
                end if.
            end loop.

            compute sort_it = make(m,m,0).

            loop i1 = 1 to m.
                loop i2 = 1 to m.
                    do if i1 <> i2.
                        do if clus_id(items(i1,2)) = clus_id(items(i2,2)) and items(i1,1) = items(i2,1).
                            compute sort_it(i1,i2) = 1.
                        end if.
                        do if clus_id(items(i1,2)) = clus_id(items(i2,2)) and items(i1,1) <> items(i2,1).
                            compute sort_it(i1,i2) = 0.
                        end if.
                        do if clus_id(items(i1,2)) <> clus_id(items(i2,2)) and items(i1,1) = items(i2,1).
                            compute sort_it(i1,i2) = 0.
                        end if.
                        do if clus_id(items(i1,2)) <> clus_id(items(i2,2)) and items(i1,1) <> items(i2,1).
                            compute sort_it(i1,i2) = 1.
                        end if.
                    end if.
                end loop.
            end loop.

            do if csum(rsum(sort_it)) = m*(m-1).
                compute sortcorr = 1.
            else if csum(rsum(sort_it)) <m*(m-1).
                compute sortcorr = 0.
            end if.

        end if.

        ******************************************
            *** COMBINE  FOR ALL SAMPLES:
        ******************************************

        do if (i_simu = 1).
            compute ovrl_all = overall.
            compute sort_all = sortcorr.
        end if.

        do if (i_simu >1).
            compute ovrl_all = {ovrl_all; overall}.
            compute sort_all = {sort_all; sortcorr}.
        end if.

    END LOOP.

    compute min_simu = cmin(ovrl_all).
    compute max_simu = cmax(ovrl_all).
    compute mn_simu = csum(ovrl_all)/n_simu.
    compute mn_help = make(n_simu, ncol(ovrl_all),0).

    loop i1 = 1 to nrow(ovrl_all).
        loop i2 = 1 to ncol(ovrl_all).
            compute mn_help(i1,i2) = mn_simu(1,i2).
        end loop.
    end loop.

    compute diff_gl = ovrl_all - mn_help.
    compute sd_simu = sqrt((cssq(diff_gl))/n_simu).

    * comparison with spacingh from data: (sim = simu, dat = data)
    * cutoff: according to z-cutoffs, 99%, one-tailed: 2.33

    compute spch_sim = mn_simu(1,1).
    compute spch_sd = sd_simu(1,1).
    compute cutoff = spch_sim + 2.33*spch_sd.
    compute comp = cutoff - spch_dat.

    compute overall = {mn_simu; sd_simu; min_simu; max_simu}.

    compute accuracy = csum(sort_all)/n_simu * 100.

    compute simu = {n_simu, n, m, p}.

    print /Title = "SAMPLES OF THE SIMULATION STUDY".

    print /Title = "The following tables show results for simulated samples from the population with perfect circumplex clusters,".

    print accuracy
        /rlabels = " Percentage"
        /clabels = "Sorting accuracy".

    print /Title = "Sorting accuracy: Percentage of simulated samples that sorted the items according to the intended circumplex clusters.".

    print overall
        /rlabels = "Mean", "SD", "Minimum", "Maximum"
        /clabels = "Spacing (with h²)", "Spacing", "Between-cluster spacing",  "Within-cluster proximity ".

    print spch_dat
        /clabels = "Spacing (with h²)"
        /rlabels = "in data".

    print /Title = "Recommendation: Circumplex fit of the empirical data is acceptable if 'spacing (with h²) in the empirical data is not larger than ".
    print /Title = "mean 'spacing (with h²)' + 2.33 SD from the simulated samples in Cluster-Circ Simu (corresponding to the cumulative".
    print /Title = "probability of the standard normal distribution for p < .01, one-tailed).".

    do if comp >= 0.
        print /Title = "Here: Empirical 'spacing (with h²)' is within mean-spc_h + 2.33 SD from Cluster-Circ-Simu --> Circumplex fit acceptable.".
    else if comp < 0.
        print /Title = "Here: Empirical 'spacing (with h²)' is larger than mean-spc_h + 2.33 SD from Cluster-Circ-Simu --> Low circumplex fit".
    end if.

     print /Title = "Range of all Cluster-Circ coefficients: 0-1 (0 = perfect circumplex spacing).".
     print  /Title = "The decision for the final item clusters is based on overall 'spacing (with h²').".
     print  /Title = "The manuscript that presents Cluster-Circ has been submitted to a peer-reviewed journal.".
     print /Title = "When using Cluster-Circ, please cite the preprint version at the current stage of the publication process:".
     print /Title = "https://psyarxiv.com/yf37w/.". 

END MATRIX.


OUTPUT SAVE NAME = *
    OUTFILE="RESULTS_CLUSTERCIRC_SIMU.spv"
    LOCK=NO.


