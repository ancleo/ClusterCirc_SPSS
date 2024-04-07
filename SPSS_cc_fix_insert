* Encoding: UTF-8.
***************************************************************************************
    *** CLUSTER-CIRC: SORTING ITEMS INTO CIRCUMPLEX-CLUSTERS
***************************************************************************************

******************************************
    *** CLUSTERCIRC-FIX: cc_fix
******************************************
    *** Computes ClusterCirc coefficients for user-defined item clusters. Items need to be sorted in the data file
    *** before analysis, e.g. cluster 1 = items 1 to 6, cluster 2 = items 7 to 10, etc.
    *** ClusterCirc-fix coefficients for user defined item clusters can be compared to ClusterCirc coefficients for item clusters found by ClusterCirc-Data.

*********************************************************************************************************************
    *** INSERTIONS: FILL IN SPECIFICATIONS OF YOUR DATA 
*********************************************************************************************************************
    *** Replace each element with "$INSERT" with the relevant information from your data
    
    *** $INSERT_path:         Working directory/path with the dataset and where the results of Cluster-Circ should be saved (line 31)
    *** $INSERT_p:               Number of clusters. Minimum = 2 (line 38)
    *** $INSERT_m:              Number of variables (line 39)
    *** $INSERT_n:               Sample size (line 40)
    *** $INSERT_data:         Name of the data file with the items, ending ".sav"; the dataset should only contain the items on which Cluster-Circ should be performed (line 50)
    *** $INSERT_limits:        Column vector with the position of the last item of each cluster, e.g {6;10;18} for 3 clusters with C1 = i1-i6, C2 = i7-i10, C3 = i11-i18 (line 109)
    *** w_com:                       Is "TRUE" if items are weighted by their communalities. If different weights should be used, change to "FALSE" (line 42).
    *** weights:                      If different weights (other than item communalities) should be used, save a data file called "weights.sav" in the working directory
    ***                                      containing a single vector with weights for the items. Weights need to be in the same order as the variables in the data and should be the same
    ***                                      weights as used in cc_data for fair comparison of spc_w (is called in line 105).
    


cd " $INSERT_path". 

new file.
dataset name DataSet1 window = front.

input program.
    loop #1 = 1 to 1.
        compute p =  $INSERT_p.
        compute m = $INSERT_m.
        compute n =  $INSERT_n.
        string w_com (a10).
        compute w_com = 'TRUE'.
        end case.
    end loop.
    end file.
end input program.
execute.
save outfile = "help.sav".

get file = "$INSERT_data".
    

********************************************************************
    *** PCA FOR INITIAL LOADING ON CIRCUMPLEX AXES
********************************************************************

OMS /SELECT ALL /DESTINATION VIEWER=NO.
FACTOR
    /VARIABLES ALL
    /MATRIX=OUT(FAC="Data_A.sav")
    /MISSING LISTWISE   /PRINT KMO  /CRITERIA FACTORS(2) ITERATE(250)
    /EXTRACTION PC  /CRITERIA ITERATE(250) NOKAISER /ROTATION NOROTATE /METHOD CORRELATION.
OMSEND.

get file = "Data_A.sav".
delete variables ROWTYPE_ FACTOR_ .
flip variables all.
delete variables CASE_LBL.
save outfile = "Data_A.sav".
get file "Data_A.sav".

output close *.

set mxloops = 300000000000.
set mdisplay TABLES.

get file = "Data_A.sav".

MATRIX.

    get A /file="Data_A.sav".
    get p /variables p /file = "help.sav".
    get m /variables m /file = "help.sav".
    get q /variables q /file = "help.sav".
    get w_com /variables w_com /file ="help.sav".

    *************************
        *** PREPARATION
    *************************
        * Kaiser-normalization of loadings for simpler computation of angles

    compute h_sq = mdiag(rssq(A)).
    compute h_rt = sqrt(h_sq).
    compute h_rtv = rsum(h_rt).
    compute hsq_mn = msum(h_sq)*(1/m).
    compute A_k = inv(h_rt)*A.

        * Define item weights:

    do if w_com = "TRUE".
        compute w = rsum(h_sq).
    end if.

    do if w_com = "FALSE".
        get w /file = "weights.sav".
    end if.

    compute w_mn = csum(w)/m.
    compute limits = $INSERT_limits.

        * Compute and sort theta: item angles in degrees
        * r = radians, p = positive
        
    compute pi=4*ARTAN(1).
    compute A_pos = abs(A_k).
    compute th_rp = arsin(A_pos(1:m,1)).
    compute th_r = make(m,1,0).
    compute theta = make(m,1,0).

        * Computation of angles depends on the quadrant (loadings positive/negative)
    
    compute item_no = make(m,1,0).

    loop i = 1 to m.
    compute item_no(i) = i.
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

        * Fixed item-cluster assignments according to 'limits':
   
    compute c_m = make(p,1,0).

    loop c = 1 to p.
        do if c = 1.
            compute c_m(c) = limits(c).
            compute ci = make(limits(c), 1, c).
            compute ci_v = ci.
        end if.
        do if c > 1.
            compute c_m(c) = limits(c) - limits(c-1).
            compute ci = make(limits(c)-limits(c-1),1,c).
            compute ci_v = {ci_v; ci}.
        end if.
    end loop.

         * Compute cluster angle as the center between the outer items in cluster

        compute c_min = make(p,1,0).
        compute c_max = make(p,1,0).
         compute c_rng = make(p,1,0).
        compute c_ang = make(p,1,0).

        loop c = 1 to p.
            compute c_min(c) = 361.
            compute c_max(c) = 0.

            loop i = 1 to m.
                do if ci_v(i) = c and theta(i) <= c_min(c).
                    compute c_min(c) = theta(i).
                end if.
                do if ci_v(i) = c and theta(i) >= c_max(c).
                    compute c_max(c) = theta(i).
                end if.
            end loop.

            do if c_max(c)-c_min(c) < 180.
                compute c_ang(c) = (c_max(c)+c_min(c))/2.
                compute c_rng(c) = c_max(c)-c_min(c).
            end if.
      
            * Special case: Cluster at approx. 0° could have a range of > 180.
            * Change item angles with help objects.

            do if c_max(c)-c_min(c) > 180.

                compute theta_h = theta.
                compute c_minh = 361.
                compute c_maxh = 0.

                loop i = 1 to m.
                    do if theta(i) > 180.
                        compute theta_h(i) = theta(i)-360.
                    end if.
                    do if ci_v(i) = c and theta_h(i) <= c_minh.
                        compute c_minh = theta_h(i).
                    end if.
                    do if ci_v(i) = c and theta_h(i) >= c_maxh.
                        compute c_maxh = theta_h(i).
                    end if.
                end loop.

                compute c_max(c) = c_maxh.
                compute c_min(c) = 360+c_minh.
                compute c_rng(c) = c_maxh -(c_minh).
                compute c_ang(c) = (c_minh + c_maxh)/2.

            end if.

        end loop.

       * Clusters and items need to be sorted before computing spacing indices

    compute ival_h = {item_no, theta, w}.
    compute cval_h = make(p,4,0).    
    
    compute c_rnk = grade(c_ang).

        loop c1 = 1 to p.
            loop c2 = 1 to p.
                do if c_rnk(c2) = c1.
                    compute cval_h(c1,1) = c2.
                    compute cval_h(c1,2) = c_m(c2).
                    compute cval_h(c1,3) = c_ang(c2).
                    compute cval_h(c1,4) = c_rng(c2).
                end if.
            end loop.
        end loop.

        compute c_no = make(m,1,0).
        compute cval = cval_h.

        loop i = 1 to m.
            loop c  = 1 to p.
                do if ci_v(i) = cval_h(c,1).
                    compute c_no(i) = c.
                    compute cval(c,1) = c.
                end if.
            end loop.
        end loop.

        * Help vectors for computation of distances:
        * Allow for negative angles (n) in first cluster
        
        compute ival_n = {c_no, ival_h}.

        loop i = 1 to m.
            do if c_no(i) = 1 and ival_h(i,2) > 180.
                compute ival_n(i,3) = ival_n(i,3) - 360.
            end if.
        end loop.

        compute rk_iang = grade(ival_n(:,3)).
        compute ival = make(m,4,0).

        loop i1 = 1 to m.
            loop i2 = 1 to m.
                do if rk_iang(i2) = i1.
                    compute ival(i1,:) = ival_n(i2,:).
                end if.
            end loop.
        end loop.

*********************************
 *** CLUSTERCIRC INDICES
*********************************

        compute ic_dis = make(m,p,0).
        compute space = 360/p.
        compute ic_dev = make(m,p,0).
        compute ic_devp = make(m,p,0).
        compute ic_dw = make(m,p,0).

        loop i = 1 to m.
            loop c1 = 1 to p.
                compute c2 = ival(i,1).
                compute i_ang =  ival(i,3).
                compute c_ang = cval(c1,3).
                compute i_w = ival(i,4).
                compute ic_dis(i,c1) = i_ang - c_ang.
                compute id_dis = (c2-1)*space - (c1-1)*space.
                compute ic_dev(i,c1) = ic_dis(i,c1) - id_dis.
                compute ic_devp(i,c1) = ic_dev(i,c1)/space.
                compute ic_dw(i,c1) = ic_devp(i,c1)* sqrt(i_w).
            end loop.
        end loop.

        * Item-spacing:

        compute ispc_sq = rssq(ic_devp)/p.
        compute ispc = sqrt(ispc_sq).

        * With weights  for spacing index (not interpretable on item level)

        compute ispc_wsq = rssq(ic_dw)/(p*w_mn).
        compute ispc_w = sqrt(ispc_wsq).

        * Overall spacing

        compute spc_sq  = csum(ispc_sq)/m.
        compute spc = sqrt(spc_sq).

        compute spc_wsq = csum(ispc_wsq)/m.
        compute spc_w = sqrt(spc_wsq).
    
        * Between-cluster spacing

        compute c_dis = make(p,p,0).
        compute space = 360/p.
        compute c_dev = make(p,p,0).
        compute c_devp = make(p,p,0).

        loop c1 = 1 to p.
            loop c2 = 1 to p.
                compute c_dis(c1,c2) = cval(c1,3)-cval(c2,3).
                compute id_dis = (c1-1)*space - (c2-1)*space.
                compute c_dev(c1,c2) = (c_dis(c1,c2) - id_dis).
                compute c_devp(c1,c2) = c_dev(c1,c2)/space.
            end loop.
        end loop.

        compute cbcs_psq = rssq(c_devp)/(p-1).
        compute cbcs_p = sqrt(cbcs_psq).

        * Overall between-cluster spacing

        compute bcs_p = sqrt(csum(cbcs_psq)/p).

        * Within-cluster proximity

        compute ic_disw = make(m,1,0).
        compute ic_disp = make(m,1,0).
        compute cwcp_sqp = make(p,1,0).
        compute cwcp_p = make(p,1,0).

        loop c = 1 to p.

            compute c_m = cval(c,2).
            compute c_icdisp = make(m,1,0).
            loop i = 1 to m.
                do if ival(i,1) = c.
                    compute ic_disw(i) = abs(ic_dis(i,c)).
                    compute ic_disp(i) = (ic_disw(i)/(360/(2*p))).
                    compute c_icdisp(i) = ic_disp(i).
                end if.
            end loop.
            compute cwcp_sqp(c) = cssq(c_icdisp)/c_m.
            compute cwcp_p(c) = sqrt(cwcp_sqp(c)).

        end loop.

        * Overall within cluster proximity

        compute wcp_p = sqrt(csum(cwcp_sqp)/p).

        * Final values

        compute items ={ival, ispc, ic_disp}.    
        compute clusters= {cval, cbcs_p, cwcp_p}.
        
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
        
        compute overall = {spc_w, spc,  bcs_p, wcp_p}.

         print  /Title = "RESULTS CLUSTER-CIRC FIX:".

        print overall
            /rlabels = " "
            /clabels = "Spacing (weighted)", "Spacing", "Between-cluster spacing",  "Within-cluster proximity ".

        print clusters
            /clabels = "Cluster", "Items", "Angle", "Range",  "Between-cluster spacing ",  "Within-cluster proximity".

        print items
            /clabels = "Cluster", "Item", "Angle", "Weight (default = communality",  "Item-cluster spacing ", "Distance to cluster center (0-1)".

        print /Title = "Range of all Cluster-Circ coefficients: 0-1 (0 = perfect circumplex spacing).".
        print  /Title = "Weights are item communalities if not otherweise specified.".
        print /Title = " ".
        print  /Title = "The manuscript that presents Cluster-Circ has been submitted to a peer-reviewed journal.".
        print /Title = "When using Cluster-Circ, please cite the preprint version at the current stage of the publication process:".
        print /Title = "(Note to reviewers: Please don't open the link below to ensure double-blind peer-review)".
        print /Title = "https://psyarxiv.com/yf37w/.". 

END MATRIX.

