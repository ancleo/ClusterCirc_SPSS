* Encoding: UTF-8.
* Encoding: .
***************************************************************************************
    *** CLUSTER-CIRC: SORTING ITEMS INTO CIRCUMPLEX-CLUSTERS
***************************************************************************************

******************************************
    *** CLUSTER-CIRC DATA: cc_data
******************************************
    *** Finds item clusters with optimal circumplex spacing in your data. 
    *** PCA without rotation will be performed before running Cluster-Circ (option in R: analyze loadings directly).
    *** cc_data can read an SPSS data-sheet (.sav) with the items that should be included into the analysis.
    *** cc_data can be performed on any number of items that is compatible with the SPSS requirements. Listwise deletion of missings if not otherwise specified.

*********************************************************************************************************************
    *** INSERTIONS: FILL IN SPECIFICATIONS OF YOUR DATA / DEMONSTRATION ON EXEMPLARY DATA
*********************************************************************************************************************
    *** Replace each element with "$INSERT" with the relevant information from your data
    
    *** $INSERT_path:         Working directory/path with the dataset and where the results of Cluster-Circ should be saved (line 34)
    *** $INSERT_p:               Desired number of clusters. Minimum = 2 (line 41)
    *** $INSERT_m:              Number of variables (line 42)
    *** $INSERT_n:               Sample size (line 43)
    *** $INSERT_data:         Name of the data file with the items, ending ".sav"; the dataset should only contain the items on which Cluster-Circ should be performed (line 56)
    *** q:                                  Precision index for the algorithm. Precision is higher for larger values. Default = 10. 
    ***                                      q can be adjusted to increase the precision of the algorithm. Must be an integer > 0 (line 45).
    *** w_com:                       Is "1" if items are weighted by their communalities. If different weights should be used, change to "0" (line 46).
    *** weights:                      If different weights (other than item communalities) should be used, save a data file called "weights.sav" in the working directory 
    ***                                       containing a single vector with weights for the items (is called in line 111).
    *** e                                   Cluster weight (0 <= e <= 1) defining the importance of within-cluster proximity versus equal cluster spacing. Default is 1/p, weighing all
    ***                                       clusters equally. e = 0: Maximum importance of between-cluster spacing, within-cluster proximity is ignored. e = 1: Maximum importance of 
    ***                                       within-cluster proximity, between-cluster spacing is ignored. (line 44)

cd "$INSERT_path".

new file.
dataset name DataSet1 window = front.

input program.
    loop #1 = 1 to 1.
        compute p =  $INSERT_p.
        compute m = $INSERT_m.
        compute n =  $INSERT_n.
        compute e = 1/p.
        compute q = 10.
        compute w_com = 1.
        end case.
    end loop.
    end file.
end input program.
execute.
save outfile = "help.sav".

get file = "help.sav".

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

**************************
    *** CLUSTER-CIRC 
**************************

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

        * Define item weights as communalities:

    do if w_com = 1.
        compute w = rsum(h_sq).
        save w /variables weights /outfile = "weights.sav".
    end if.

END MATRIX.


set mxloops = 300000000000.
set mdisplay TABLES.

get file = "Data_A.sav".

MATRIX.

    get A /file="Data_A.sav".
    get p /variables p /file = "help.sav".
    get m /variables m /file = "help.sav".
    get q /variables q /file = "help.sav".
    get e /variables e /file = "help.sav".
    get w /variables weights /file ="weights.sav".

*************************
        *** PREPARATION
    *************************
        * Kaiser-normalization of loadings for simpler computation of angles

    compute h_sq = mdiag(rssq(A)).
    compute h_rt = sqrt(h_sq).
    compute h_rtv = rsum(h_rt).
    compute hsq_mn = msum(h_sq)*(1/m).
    compute A_k = inv(h_rt)*A.

    compute w_mn = csum(w)/m.

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
                compute ival(i2,3) = w(i1).
            end if.
        end loop.
    end loop.

    ***********************************************
        *** CLUSTERCIRC ALGORITHM
    ***********************************************
    
    compute spacingw = 361.

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

        * Compute spacing_w for each division

        compute ic_dis = make(m,p,0).
        compute space = 360/p.
        compute ic_dev = make(m,p,0).
        compute ic_devp = make(m,p,0).
        compute ic_dw = make(m,p,0).
        compute ic_dwe = make(m,p,0).

        loop i = 1 to m.
            loop c1 = 1 to p.

                compute c2 = ival_h(i,1).
                compute i_ang =  ival_h(i,3).
                compute c_ang = cvalh(c1,3).
                compute i_w = ival_h(i,4).
                compute e_own = e.
                compute e_others = (1-e)/(p-1).
                compute ic_dis(i,c1) = i_ang - c_ang.
                compute id_dis = (c2-1)*space - (c1-1)*space.
                compute ic_dev(i,c1) = ic_dis(i,c1) - id_dis.
                compute ic_devp(i,c1) = ic_dev(i,c1)/space.
                compute ic_dw(i,c1) = ic_devp(i,c1)* sqrt(i_w).

                do if c1 = c2.
                    compute ic_dwe(i,c1) = ic_dw(i,c1)*sqrt(e_own).
                end if.

                do if c1 <> c2.
                    compute ic_dwe(i,c1) = ic_dw(i,c1)*sqrt(e_others).
                end if.

            end loop.
        end loop.

        * Item-spacing:

        compute ispc_sq = rssq(ic_devp)/p.
        compute ispc = sqrt(ispc_sq).

        * With weights for spacing index (not interpretable on item level)

        compute ispc_wsq = rssq(ic_dwe)/ w_mn.
        compute ispc_w = sqrt(ispc_wsq).

        * Overall spacing

        compute spc_sq  = csum(ispc_sq)/m.
        compute spc = sqrt(spc_sq).

        compute spc_wsq = csum(ispc_wsq)/m.
        compute spc_w = sqrt(spc_wsq).

        * Dismiss partitions with empty clusters by making spc_w larger than
        * the initial spacingw (361). If this happens for all possible divisions,
        * the number of clusters is too large.

        loop c = 1 to p.
            do if c_m(c) = 99.
                compute spc_w = 50000.
            end if.
        end loop.

        do if spc_w < spacingw.
            compute spacingw = spc_w.
            compute spacing  = spc.
            compute items = {ival_h, ispc}.
            compute clusters = cvalh.
            compute ic_dist = ic_dis.
        end if.

    END LOOP.

    do if spacingw = 361.
        print/Title "Cluster-Circ could not finish, at least one of the clusters is empty. Try a smaller number of clusters or include more variables.".

    else if spacingw < 361.

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
        
        compute overall = {spacingw, spacing,  bcs_p, wcp_p}.
    
        * Parameters for cc_simu

        compute c_wrange = clusters(:,4).
        compute wrange_c = csum(c_wrange)/p.

        save {hsq_mn, wrange_c, spacingw, e}  /outfile = "help_CC_data.sav"
            /variables hsq_mn wrange_c spacingw e.

        print  /Title = "RESULTS CLUSTER-CIRC DATA:".

        print overall
            /rlabels = " "
            /clabels = "Spacing (weighted)", "Spacing", "Between-cluster spacing",  "Within-cluster proximity ".

        print clusters
            /clabels = "Cluster", "Items", "Angle", "Range",  "Between-cluster spacing ",  "Within-cluster proximity".

        print items
            /clabels = "Cluster", "Item", "Angle", "Weight (default = communality",  "Item-cluster spacing ", "Distance to cluster center (0-1)".

        print /Title = "Range of all Cluster-Circ coefficients: 0-1 (0 = perfect circumplex spacing).".
        print  /Title = "The decision for the final item clusters is based on overall 'spacing (weighted)'.".
        print  /Title = "Weights are item communalities if not otherweise specified to account for (un-)reliability of measurement.".
        print  /Title = "The manuscript that presents Cluster-Circ has been submitted to a peer-reviewed journal.".
        print /Title = "When using Cluster-Circ, please cite the preprint version at the current stage of the publication process:".
        print /Title = "(Note to reviewers: Please don't open the link below to ensure double-blind peer-review)".
        print /Title = "https://psyarxiv.com/yf37w/.". 

    end if.

END MATRIX.

OUTPUT SAVE NAME = *
    OUTFILE="RESULTS_CLUSTERCIRC.spv"
    LOCK=NO.
