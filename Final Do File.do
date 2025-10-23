***************************************************************************
**************ALSPAC BODY COMPOSITION AND CV PHENOTYPES********************
*************************Authors: Mok & Iyer*******************************
*********************************2025**************************************
***************************************************************************

version 19
clear all
macro drop _all
set more off
set maxvar 10000


//////////////////////////////////////////////////////////////////////////
/////////////////////////DATA CLEANING AND PREP///////////////////////////
//////////////////////////////////////////////////////////////////////////


cd "S:\ICS_Student_Projects\2021-22\AiL Paper\Paper\Analysis"
use "ALSPAC+HR merged.dta", clear
drop _merge
merge m:1 cidB3076 qlet using "S:\ICS_Student_Projects\2021-22\AiL Paper\archive\ALSPAC+HR merged Jul72023.dta", keepusing(brshr)

duplicates drop

//Covariates//

clonevar sex=kz021
replace sex=. if kz021<0

clonevar ses=c755
replace ses=. if c755<0

tab ses
tab ses, nolabel
replace ses = . if ses < 0
replace ses = . if ses == 65 //remove class "Armed forces//

clonevar physact=fh5015
replace physact=. if fh5015<0 

clonevar geneticscore = SCORE

//Tissue Type Exposures//

gen fm=tfvoxel*(0.918/1000)
gen ffm=wt-fm
gen ffmi=ffm/ht^2
gen fmi=fm/ht^2
gen bf=(fm/wt)*100

//Cardiac Structure Outcomes//

clonevar lvm=mrilvmass
gen lvmiffm=lvm/ffm

//Cardiac Function Outcomes//

clonevar edv=mrilvedv
clonevar esv=mrilvesv
clonevar sv=mrilvsv
clonevar co=mriao1co
clonevar hr=brshr
clonevar ef=mrilvef
gen cmw = mriao1co * ((dbppwa*2 + pwaasp)/3)
gen lvws=pwaasp*(1+(3*((esv/ht^2.7)/(lvm/ht^2.7))))

//Blood Pressure Outcomes//

clonevar psbp = sbppwa //peripheral sbp//
clonevar dbp = dbppwa
generate ppp = sbppwa-dbppwa //peripheral pulse pressure//
clonevar csbp = pwaasp //central sbp//
generate cpp = csbp - dbp
generate map = (dbppwa*2 + sbppwa)/3 //peripheral mean arterial pressure//

//Autonomic Function Outcomes//

clonevar hrvsdnn = brstdhrv //heart rate variability sdnn//
clonevar hrvrmsdd = brstdhrv2 //heart rate variability msdd//
clonevar hrvti = brstdhrv3 //heart rate variability ti//
clonevar brs = brsslope //baroreflex sensitivity//
clonevar svr = mriao1svr

//Drop Four Physiologically Impossible/Improbable Outliers//

drop if bf > 70
drop if cidB3076 == 15927
drop if cidB3076 == 8695
drop if cidB3076 == 10243

//Drop Those Missing Body Composition Measures//

keep if fm!=.

//Order Variables//

order bmi bf wt ht fm ffm ffmi fmi lvm lvmiffm lvws psbp ppp csbp dbp cpp map edv esv ef sv co hrvsdnn hrvrmsdd hrvti brs hr svr age sex ses physact geneticscore 
zscore bmi bf wt ht fm ffm ffmi fmi lvm lvmiffm lvws psbp ppp csbp dbp cpp map edv esv ef sv co hrvsdnn hrvrmsdd hrvti brs hr svr age physact 

save "Data_for_analysis.dta", replace

//////////////////////////////////////////////////////////////////////////
///////////////////////////ANALYSIS///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


*******************CHARACTERISTICS FOR TABLE 1****************************

//Exposures//

sum age bmi bf fm ffm, det
tab sex

//Outcomes//

sum lvm lvws 
sum csbp cpp dbp map
sum edv esv sv ef co cmw
sum hrvsdnn hrvrmsdd hrvti brs svr brshr

//Covariates//

tab ses
sum physact 
tab geneticscore 

//Correlations Between Outcomes To Justify Lack of Multiple Comparisons in Analyses//
//Heatmap for these created later in supplementary section//

pwcorr lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs, sig


******************SCATTER PLOTS FOR FIGURE 1*****************************

twoway ///
    (scatter fm ffm if sex == 1, mcolor("0 0 153") msymbol(circle)) ///
    (scatter fm ffm if sex == 2, mcolor("153 0 0") msymbol(circle)), ///
    legend(order(1 "Male" 2 "Female")) xtitle("FFM (kg)") ytitle("FM (kg)")
graph export Figure1L.png, replace

twoway ///
    (scatter bf bmi if sex == 1, mcolor("0 0 153") msymbol(circle)) ///
    (scatter bf bmi if sex == 2, mcolor("153 0 0") msymbol(circle)), ///
    legend(order(1 "Male" 2 "Female")) xtitle("BMI (kg/m2)") ytitle("Body Fat (%)")
graph export Figure1R.png, replace


///////////////////////////////////////////////////////////////////
//////////////////////IMPUTING DATA////////////////////////////////
///////////////////////////////////////////////////////////////////

use "Data_for_analysis.dta", clear

recode geneticscore (1=0) (2=1)

mi set wide
mi register regular z_fm z_ffm sex age 
mi register imputed lvm edv esv sv co ef brshr lvws cmw csbp dbp map svr hrvsdnn hrvrmsdd hrvti brs ses physact geneticscore
mi update
mi impute chained (regress) lvm edv esv sv co ef brshr lvws cmw csbp dbp map svr hrvsdnn hrvrmsdd hrvti brs physact (ologit) ses (logit) geneticscore = z_fm z_ffm age i.sex, add(100) rseed(54321) dots

save "Imputed_data_for_analysis.dta", replace

*******************MI REGRESSIONS FOR TABLE 2****************************

//Below code runs four models with increasing levels of adjustment on all outcomes and exports each individually to a .csv//
//These are then modified with the appropriate labels, saved as data files, and then all appended together to form the final table//

preserve

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
	
eststo clear
eststo m1_fm_`out': quietly mi estimate, post: regress `out' z_fm 
eststo m2_fm_`out': quietly mi estimate, post: regress `out' z_fm age i.sex
eststo m3_fm_`out': quietly mi estimate, post: regress `out' z_fm age i.sex physact i.ses i.geneticscore
eststo m4_fm_`out': quietly mi estimate, post: regress `out' z_fm age i.sex physact i.ses i.geneticscore z_ffm

esttab m1_fm_`out' m2_fm_`out' m3_fm_`out' m4_fm_`out' ///
       using `out'_fm_results_imputed.csv, replace ///
       cells("b(fmt(%9.2f)) ci_l(fmt(%9.2f)) ci_u(fmt(%9.2f)) p(fmt(3))") ///
       label nodepvar ///
       varlabels(fm "FM") ///
	   keep(z_fm) ///
       nonumber alignment(D) noobs compress ///
	   plain
	   	   
}

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
	
eststo clear
eststo m1_ffm_`out': quietly mi estimate, post: regress `out' z_ffm 
eststo m2_ffm_`out': quietly mi estimate, post: regress `out' z_ffm age i.sex
eststo m3_ffm_`out': quietly mi estimate, post: regress `out' z_ffm age i.sex physact i.ses i.geneticscore
eststo m4_ffm_`out': quietly mi estimate, post: regress `out' z_ffm age i.sex physact i.ses i.geneticscore z_fm

esttab m1_ffm_`out' m2_ffm_`out' m3_ffm_`out' m4_ffm_`out' ///
       using `out'_ffm_results_imputed.csv, replace ///
       cells("b(fmt(%9.2f)) ci_l(fmt(%9.2f)) ci_u(fmt(%9.2f)) p(fmt(3))") ///
       label nodepvar ///
       varlabels(ffm "FFM") ///
	   keep(z_ffm) ///
       nonumber alignment(D) noobs compress ///
	   plain
	   
}

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
import delimited "`out'_fm_results_imputed.csv", rowrange(3:3) colrange(2:) clear
rename m1_fm_`out' m1_beta
rename v2 m1_lci
rename v3 m1_uci
rename v4 m1_p
rename m2_fm_`out' m2_beta
rename v6 m2_lci
rename v7 m2_uci
rename v8 m2_p
rename m3_fm_`out' m3_beta
rename v10 m3_lci
rename v11 m3_uci
rename v12 m3_p
rename m4_fm_`out' m4_beta
rename v14 m4_lci
rename v15 m4_uci
rename v16 m4_p
save "`out'_fm_results_imputed.dta", replace
}

use "Data_for_analysis.dta", clear

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
import delimited "`out'_ffm_results_imputed.csv", rowrange(3:3) colrange(2:) clear
rename m1_ffm_`out' m1_beta
rename v2 m1_lci
rename v3 m1_uci
rename v4 m1_p
rename m2_ffm_`out' m2_beta
rename v6 m2_lci
rename v7 m2_uci
rename v8 m2_p
rename m3_ffm_`out' m3_beta
rename v10 m3_lci
rename v11 m3_uci
rename v12 m3_p
rename m4_ffm_`out' m4_beta
rename v14 m4_lci
rename v15 m4_uci
rename v16 m4_p
save "`out'_ffm_results_imputed.dta", replace
}


use "lvm_fm_results_imputed.dta", clear
append using "lvm_ffm_results_imputed.dta", force
append using "edv_fm_results_imputed.dta", force
append using "edv_ffm_results_imputed.dta", force
append using "esv_fm_results_imputed.dta", force
append using "esv_ffm_results_imputed.dta", force
append using "sv_fm_results_imputed.dta", force
append using "sv_ffm_results_imputed.dta", force
append using "ef_fm_results_imputed.dta", force
append using "ef_ffm_results_imputed.dta", force
append using "co_fm_results_imputed.dta", force
append using "co_ffm_results_imputed.dta", force
append using "hr_fm_results_imputed.dta", force
append using "hr_ffm_results_imputed.dta", force
append using "lvws_fm_results_imputed.dta", force
append using "lvws_ffm_results_imputed.dta", force
append using "cmw_fm_results_imputed.dta", force
append using "cmw_ffm_results_imputed.dta", force
append using "csbp_fm_results_imputed.dta", force
append using "csbp_ffm_results_imputed.dta", force
append using "cpp_fm_results_imputed.dta", force
append using "cpp_ffm_results_imputed.dta", force
append using "dbp_fm_results_imputed.dta", force
append using "dbp_ffm_results_imputed.dta", force
append using "map_fm_results_imputed.dta", force
append using "map_ffm_results_imputed.dta", force
append using "svr_fm_results_imputed.dta", force
append using "svr_ffm_results_imputed.dta", force
append using "hrvsdnn_fm_results_imputed.dta", force
append using "hrvsdnn_ffm_results_imputed.dta", force
append using "hrvrmsdd_fm_results_imputed.dta", force
append using "hrvrmsdd_ffm_results_imputed.dta", force
append using "hrvti_fm_results_imputed.dta", force
append using "hrvti_ffm_results_imputed.dta", force
append using "brs_fm_results_imputed.dta", force
append using "brs_ffm_results_imputed.dta", force

save "Table_2_Results.dta", replace
export delimited "Table_2_Results.csv", replace

restore

************************GRAPHS FOR FIGURE 3*****************************

foreach out of varlist(co sbp cmw lvm sv hr) {

quietly mi estimate: regress `out' c.bmi##c.bf i.sex age i.ses physact i.geneticscore 
quietly mimrgns, dydx(bf) at(bmi=(25(5)35)) vsquish cmdmargins
quietly mimrgns, at(bf=(15 25 35) bmi=(25(5)35)) vsquish cmdmargins
quietly marginsplot, noci x(bmi) recast(line) xlabel(25(5)35) xtitle("BMI (kg/m2)") ytitle("Left Ventricular Mass (g)") legend(off) ///
    plot1opts(lcolor("153 0 0")) ///
    plot2opts(lcolor("0 0 153")) ///
    plot3opts(lcolor("0 102 0")) ///
graph export "`out'_graph.png", replace 

}


/////////////////////////////////////////////////////////////////////////
//////////////////////SUPPLEMENTARY ANALYSES/////////////////////////////
/////////////////////////////////////////////////////////////////////////

use "Data_for_analysis.dta", clear


******************CORRELATION HEATMAP FOR OUTCOMES***********************

correlate lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs
matrix r = r(C)

heatplot r, color(hcl diverging, intensity(0.6)) ///
    xlabel(, angle(45) labsize(small)) ///
    ylabel(, labsize(small)) ///
	keylabels(, format(%4.1f)) ///
    aspectratio(1) ///
    legend(label(1 "-1.00") label(2 "-0.75") label(3 "-0.50") label(4 "-0.25") label(5 "0.00") label(6 "0.25") label(7 "0.50") label(8 "0.75") label(9 "1.00")) ///
    xsize(20) ysize(20) ///
	lower
graph export "Heatmap.png", replace

***********************DATA FOR SUPP TABLE 1****************************

misstable sum z_bmi z_bf z_fm z_ffm sex age lvm edv esv sv co ef brshr lvws cmw csbp dbp map svr hrvsdnn hrvrmsdd hrvti brs ses physact geneticscore


***********************COMPLETE CASE ANALYSES****************************

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
	
eststo clear
eststo m1_fm_`out': quietly regress `out' z_fm 
eststo m2_fm_`out': quietly regress `out' z_fm age i.sex
eststo m3_fm_`out': quietly regress `out' z_fm age i.sex physact i.ses i.geneticscore
eststo m4_fm_`out': quietly regress `out' z_fm age i.sex physact i.ses i.geneticscore z_ffm

esttab m1_fm_`out' m2_fm_`out' m3_fm_`out' m4_fm_`out' ///
       using `out'_fm_results.csv, replace ///
       cells("b(fmt(%9.2f)) ci_l(fmt(%9.2f)) ci_u(fmt(%9.2f)) p(fmt(3))") ///
       label nodepvar ///
       varlabels(fm "FM") ///
	   keep(z_fm) ///
       nonumber alignment(D) noobs compress ///
	   plain
	   	   
}

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
	
eststo clear
eststo m1_ffm_`out': quietly regress `out' z_ffm 
eststo m2_ffm_`out': quietly regress `out' z_ffm age i.sex
eststo m3_ffm_`out': quietly regress `out' z_ffm age i.sex physact i.ses i.geneticscore
eststo m4_ffm_`out': quietly regress `out' z_ffm age i.sex physact i.ses i.geneticscore z_fm

esttab m1_ffm_`out' m2_ffm_`out' m3_ffm_`out' m4_ffm_`out' ///
       using `out'_ffm_results.csv, replace ///
       cells("b(fmt(%9.2f)) ci_l(fmt(%9.2f)) ci_u(fmt(%9.2f)) p(fmt(3))") ///
       label nodepvar ///
       varlabels(ffm "FFM") ///
	   keep(z_ffm) ///
       nonumber alignment(D) noobs compress ///
	   plain
	   
}

//Create Combined Table//

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
import delimited "`out'_fm_results.csv", rowrange(3:3) colrange(2:) clear
rename m1_fm_`out' m1_beta
rename v2 m1_lci
rename v3 m1_uci
rename v4 m1_p
rename m2_fm_`out' m2_beta
rename v6 m2_lci
rename v7 m2_uci
rename v8 m2_p
rename m3_fm_`out' m3_beta
rename v10 m3_lci
rename v11 m3_uci
rename v12 m3_p
rename m4_fm_`out' m4_beta
rename v14 m4_lci
rename v15 m4_uci
rename v16 m4_p
save "`out'_fm_results.dta", replace
}

use "Data_for_analysis.dta", clear

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
import delimited "`out'_ffm_results.csv", rowrange(3:3) colrange(2:) clear
rename m1_ffm_`out' m1_beta
rename v2 m1_lci
rename v3 m1_uci
rename v4 m1_p
rename m2_ffm_`out' m2_beta
rename v6 m2_lci
rename v7 m2_uci
rename v8 m2_p
rename m3_ffm_`out' m3_beta
rename v10 m3_lci
rename v11 m3_uci
rename v12 m3_p
rename m4_ffm_`out' m4_beta
rename v14 m4_lci
rename v15 m4_uci
rename v16 m4_p
save "`out'_ffm_results.dta", replace
}

use "lvm_fm_results.dta", clear
append using "lvm_ffm_results.dta", force
append using "edv_fm_results.dta", force
append using "edv_ffm_results.dta", force
append using "esv_fm_results.dta", force
append using "esv_ffm_results.dta", force
append using "sv_fm_results.dta", force
append using "sv_ffm_results.dta", force
append using "ef_fm_results.dta", force
append using "ef_ffm_results.dta", force
append using "co_fm_results.dta", force
append using "co_ffm_results.dta", force
append using "hr_fm_results.dta", force
append using "hr_ffm_results.dta", force
append using "lvws_fm_results.dta", force
append using "lvws_ffm_results.dta", force
append using "cmw_fm_results.dta", force
append using "cmw_ffm_results.dta", force
append using "csbp_fm_results.dta", force
append using "csbp_ffm_results.dta", force
append using "cpp_fm_results.dta", force
append using "cpp_ffm_results.dta", force
append using "dbp_fm_results.dta", force
append using "dbp_ffm_results.dta", force
append using "map_fm_results.dta", force
append using "map_ffm_results.dta", force
append using "svr_fm_results.dta", force
append using "svr_ffm_results.dta", force
append using "hrvsdnn_fm_results.dta", force
append using "hrvsdnn_ffm_results.dta", force
append using "hrvrmsdd_fm_results.dta", force
append using "hrvrmsdd_ffm_results.dta", force
append using "hrvti_fm_results.dta", force
append using "hrvti_ffm_results.dta", force
append using "brs_fm_results.dta", force
append using "brs_ffm_results.dta", force

save "Supp_Table_1_Results.dta", replace
export delimited "Supp_Table_1_Results.csv", replace


*********************SEX INTERACTIONS*********************

use "Imputed_data_for_analysis.dta", replace

preserve 

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
	
eststo clear
eststo m_fm_`out': quietly mi estimate, post: regress `out' z_fm age physact i.ses i.geneticscore z_ffm if sex == 1
eststo f_fm_`out': quietly mi estimate, post: regress `out' z_fm age physact i.ses i.geneticscore z_ffm if sex == 2

esttab m_fm_`out' f_fm_`out' ///
       using `out'_fm_sex_results_imputed.csv, replace ///
       cells("b(fmt(%9.2f)) ci_l(fmt(%9.2f)) ci_u(fmt(%9.2f)) p(fmt(3))") ///
       label nodepvar ///
       varlabels(fm "FM") ///
	   keep(z_fm) ///
       nonumber alignment(D) noobs compress ///
	   plain
	   	   
}

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
	
eststo clear
eststo m_ffm_`out': quietly mi estimate, post: regress `out' z_ffm age physact i.ses i.geneticscore z_fm if sex == 1
eststo f_ffm_`out': quietly mi estimate, post: regress `out' z_ffm age physact i.ses i.geneticscore z_fm if sex == 2

esttab m_ffm_`out' f_ffm_`out' ///
       using `out'_ffm_sex_results_imputed.csv, replace ///
       cells("b(fmt(%9.2f)) ci_l(fmt(%9.2f)) ci_u(fmt(%9.2f)) p(fmt(3))") ///
       label nodepvar ///
       varlabels(ffm "FFM") ///
	   keep(z_ffm) ///
       nonumber alignment(D) noobs compress ///
	   plain
	   	   
}

//Create Combined Table//

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
import delimited "`out'_fm_sex_results_imputed.csv", rowrange(3:3) colrange(2:) clear
rename m_fm_`out' m_beta
rename v2 m1_lci
rename v3 m1_uci
rename v4 m1_p
rename f_fm_`out' f_beta
rename v6 m2_lci
rename v7 m2_uci
rename v8 m2_p

save "`out'_fm_sex_results_imputed.dta", replace
}

use "Imputed_data_for_analysis.dta", clear

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
import delimited "`out'_ffm_sex_results_imputed.csv", rowrange(3:3) colrange(2:) clear
rename m_ffm_`out' m_beta
rename v2 m1_lci
rename v3 m1_uci
rename v4 m1_p
rename f_ffm_`out' f_beta
rename v6 m2_lci
rename v7 m2_uci
rename v8 m2_p

save "`out'_ffm_sex_results_imputed.dta", replace
}


use "lvm_fm_sex_results_imputed.dta", clear
append using "lvm_ffm_sex_results_imputed.dta", force
append using "edv_fm_sex_results_imputed.dta", force
append using "edv_ffm_sex_results_imputed.dta", force
append using "esv_fm_sex_results_imputed.dta", force
append using "esv_ffm_sex_results_imputed.dta", force
append using "sv_fm_sex_results_imputed.dta", force
append using "sv_ffm_sex_results_imputed.dta", force
append using "ef_fm_sex_results_imputed.dta", force
append using "ef_ffm_sex_results_imputed.dta", force
append using "co_fm_sex_results_imputed.dta", force
append using "co_ffm_sex_results_imputed.dta", force
append using "hr_fm_sex_results_imputed.dta", force
append using "hr_ffm_sex_results_imputed.dta", force
append using "lvws_fm_sex_results_imputed.dta", force
append using "lvws_ffm_sex_results_imputed.dta", force
append using "cmw_fm_sex_results_imputed.dta", force
append using "cmw_ffm_sex_results_imputed.dta", force
append using "csbp_fm_sex_results_imputed.dta", force
append using "csbp_ffm_sex_results_imputed.dta", force
append using "cpp_fm_sex_results_imputed.dta", force
append using "cpp_ffm_sex_results_imputed.dta", force
append using "dbp_fm_sex_results_imputed.dta", force
append using "dbp_ffm_sex_results_imputed.dta", force
append using "map_fm_sex_results_imputed.dta", force
append using "map_ffm_sex_results_imputed.dta", force
append using "svr_fm_sex_results_imputed.dta", force
append using "svr_ffm_sex_results_imputed.dta", force
append using "hrvsdnn_fm_sex_results_imputed.dta", force
append using "hrvsdnn_ffm_sex_results_imputed.dta", force
append using "hrvrmsdd_fm_sex_results_imputed.dta", force
append using "hrvrmsdd_ffm_sex_results_imputed.dta", force
append using "hrvti_fm_sex_results_imputed.dta", force
append using "hrvti_ffm_sex_results_imputed.dta", force
append using "brs_fm_sex_results_imputed.dta", force
append using "brs_ffm_sex_results_imputed.dta", force

save "Supp_Table_3_Results.dta", replace
export delimited "Supp_Table_3_Results.csv", replace

restore

//Get p-values for Interaction tests//

preserve 

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
	
eststo clear
eststo i_fm_`out': quietly mi estimate, post: regress `out' c.z_fm##i.sex age physact i.ses i.geneticscore z_ffm

esttab i_fm_`out' ///
       using `out'_fm_sex_int_pvalues_imputed.csv, replace ///
       cells("p(fmt(3))") ///
       label nodepvar ///
       varlabels(fm "FM") ///
	   keep(2.sex#c.z_fm) ///
       nonumber alignment(D) noobs compress ///
	   plain
	   	   
}

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
	
eststo clear
eststo i_ffm_`out': quietly mi estimate, post: regress `out' c.z_ffm##i.sex age physact i.ses i.geneticscore z_fm

esttab i_ffm_`out' ///
       using `out'_ffm_sex_int_pvalues_imputed.csv, replace ///
       cells("p(fmt(3))") ///
       label nodepvar ///
       varlabels(ffm "FFM") ///
	   keep(2.sex#c.z_ffm) ///
       nonumber alignment(D) noobs compress ///
	   plain
	   	   
}

restore

//Create Combined Table//

preserve

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
import delimited "`out'_fm_sex_int_pvalues_imputed.csv", rowrange(3:3) colrange(2:) clear
rename i_fm_`out' p
save "`out'_fm_sex_int_pvalues_imputed.dta", replace
}

restore 

preserve

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
import delimited "`out'_ffm_sex_int_pvalues_imputed.csv", rowrange(3:3) colrange(2:) clear
rename i_ffm_`out' p
save "`out'_ffm_sex_int_pvalues_imputed.dta", replace
}

restore

use "lvm_fm_sex_int_pvalues_imputed.dta", clear

append using "lvm_ffm_sex_int_pvalues_imputed.dta", force
append using "edv_fm_sex_int_pvalues_imputed.dta", force
append using "edv_ffm_sex_int_pvalues_imputed.dta", force
append using "esv_fm_sex_int_pvalues_imputed.dta", force
append using "esv_ffm_sex_int_pvalues_imputed.dta", force
append using "sv_fm_sex_int_pvalues_imputed.dta", force
append using "sv_ffm_sex_int_pvalues_imputed.dta", force
append using "ef_fm_sex_int_pvalues_imputed.dta", force
append using "ef_ffm_sex_int_pvalues_imputed.dta", force
append using "co_fm_sex_int_pvalues_imputed.dta", force
append using "co_ffm_sex_int_pvalues_imputed.dta", force
append using "hr_fm_sex_int_pvalues_imputed.dta", force
append using "hr_ffm_sex_int_pvalues_imputed.dta", force
append using "lvws_fm_sex_int_pvalues_imputed.dta", force
append using "lvws_ffm_sex_int_pvalues_imputed.dta", force
append using "cmw_fm_sex_int_pvalues_imputed.dta", force
append using "cmw_ffm_sex_int_pvalues_imputed.dta", force
append using "csbp_fm_sex_int_pvalues_imputed.dta", force
append using "csbp_ffm_sex_int_pvalues_imputed.dta", force
append using "cpp_fm_sex_int_pvalues_imputed.dta", force
append using "cpp_ffm_sex_int_pvalues_imputed.dta", force
append using "dbp_fm_sex_int_pvalues_imputed.dta", force
append using "dbp_ffm_sex_int_pvalues_imputed.dta", force
append using "map_fm_sex_int_pvalues_imputed.dta", force
append using "map_ffm_sex_int_pvalues_imputed.dta", force
append using "svr_fm_sex_int_pvalues_imputed.dta", force
append using "svr_ffm_sex_int_pvalues_imputed.dta", force
append using "hrvsdnn_fm_sex_int_pvalues_imputed.dta", force
append using "hrvsdnn_ffm_sex_int_pvalues_imputed.dta", force
append using "hrvrmsdd_fm_sex_int_pvalues_imputed.dta", force
append using "hrvrmsdd_ffm_sex_int_pvalues_imputed.dta", force
append using "hrvti_fm_sex_int_pvalues_imputed.dta", force
append using "hrvti_ffm_sex_int_pvalues_imputed.dta", force
append using "brs_fm_sex_int_pvalues_imputed.dta", force
append using "brs_ffm_sex_int_pvalues_imputed.dta", force

save "Supp_Table_3_int_pvalues.dta", replace
export delimited "Supp_Table_3_int_pvalues.csv", replace


**********************GENETIC INTERACTIONS**************************

use "Data_for_analysis.dta", clear

mi set wide
mi register regular z_fm z_ffm sex age
mi register imputed lvm edv esv sv co ef brshr lvws cmw csbp dbp map svr hrvsdnn hrvrmsdd hrvti brs ses physact 
mi update
mi impute chained (regress) lvm edv esv sv co ef brshr lvws cmw csbp dbp map svr hrvsdnn hrvrmsdd hrvti brs physact (ologit) ses = z_fm z_ffm age i.sex, add(100) rseed(54321) dots

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
	
eststo clear
eststo m_fm_`out': quietly mi estimate, post: regress `out' z_fm age physact i.ses i.sex z_ffm if geneticscore == 1
eststo f_fm_`out': quietly mi estimate, post: regress `out' z_fm age physact i.ses i.sex z_ffm if geneticscore == 2

esttab m_fm_`out' f_fm_`out' ///
       using `out'_fm_gene_results_imputed.csv, replace ///
       cells("b(fmt(%9.2f)) ci_l(fmt(%9.2f)) ci_u(fmt(%9.2f)) p(fmt(3))") ///
       label nodepvar ///
       varlabels(fm "FM") ///
	   keep(z_fm) ///
       nonumber alignment(D) noobs compress ///
	   plain
	   	   
}

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
	
eststo clear
eststo m_ffm_`out': quietly mi estimate, post: regress `out' z_ffm age physact i.ses i.sex z_fm if geneticscore == 1
eststo f_ffm_`out': quietly mi estimate, post: regress `out' z_ffm age physact i.ses i.sex z_fm if geneticscore == 2

esttab m_ffm_`out' f_ffm_`out' ///
       using `out'_ffm_gene_results_imputed.csv, replace ///
       cells("b(fmt(%9.2f)) ci_l(fmt(%9.2f)) ci_u(fmt(%9.2f)) p(fmt(3))") ///
       label nodepvar ///
       varlabels(ffm "FFM") ///
	   keep(z_ffm) ///
       nonumber alignment(D) noobs compress ///
	   plain
	   	   
}

//Create Combined Table//

preserve

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
import delimited "`out'_fm_gene_results_imputed.csv", rowrange(3:3) colrange(2:) clear
rename m_fm_`out' l_beta
rename v2 m1_lci
rename v3 m1_uci
rename v4 m1_p
rename f_fm_`out' h_beta
rename v6 m2_lci
rename v7 m2_uci
rename v8 m2_p

save "`out'_fm_gene_results_imputed.dta", replace
}

restore

preserve

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
import delimited "`out'_ffm_gene_results_imputed.csv", rowrange(3:3) colrange(2:) clear
rename m_ffm_`out' l_beta
rename v2 m1_lci
rename v3 m1_uci
rename v4 m1_p
rename f_ffm_`out' h_beta
rename v6 m2_lci
rename v7 m2_uci
rename v8 m2_p

save "`out'_ffm_gene_results_imputed.dta", replace
}

restore

preserve

use "lvm_fm_gene_results_imputed.dta", clear
append using "lvm_ffm_gene_results_imputed.dta", force
append using "edv_fm_gene_results_imputed.dta", force
append using "edv_ffm_gene_results_imputed.dta", force
append using "esv_fm_gene_results_imputed.dta", force
append using "esv_ffm_gene_results_imputed.dta", force
append using "sv_fm_gene_results_imputed.dta", force
append using "sv_ffm_gene_results_imputed.dta", force
append using "ef_fm_gene_results_imputed.dta", force
append using "ef_ffm_gene_results_imputed.dta", force
append using "co_fm_gene_results_imputed.dta", force
append using "co_ffm_gene_results_imputed.dta", force
append using "hr_fm_gene_results_imputed.dta", force
append using "hr_ffm_gene_results_imputed.dta", force
append using "lvws_fm_gene_results_imputed.dta", force
append using "lvws_ffm_gene_results_imputed.dta", force
append using "cmw_fm_gene_results_imputed.dta", force
append using "cmw_ffm_gene_results_imputed.dta", force
append using "csbp_fm_gene_results_imputed.dta", force
append using "csbp_ffm_gene_results_imputed.dta", force
append using "cpp_fm_gene_results_imputed.dta", force
append using "cpp_ffm_gene_results_imputed.dta", force
append using "dbp_fm_gene_results_imputed.dta", force
append using "dbp_ffm_gene_results_imputed.dta", force
append using "map_fm_gene_results_imputed.dta", force
append using "map_ffm_gene_results_imputed.dta", force
append using "svr_fm_gene_results_imputed.dta", force
append using "svr_ffm_gene_results_imputed.dta", force
append using "hrvsdnn_fm_gene_results_imputed.dta", force
append using "hrvsdnn_ffm_gene_results_imputed.dta", force
append using "hrvrmsdd_fm_gene_results_imputed.dta", force
append using "hrvrmsdd_ffm_gene_results_imputed.dta", force
append using "hrvti_fm_gene_results_imputed.dta", force
append using "hrvti_ffm_gene_results_imputed.dta", force
append using "brs_fm_gene_results_imputed.dta", force
append using "brs_ffm_gene_results_imputed.dta", force

save "Supp_Table_4_Results.dta", replace
export delimited "Supp_Table_4_Results.csv", replace

restore

//Get p-values for Interaction tests//

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
	
eststo clear
eststo i_fm_`out': quietly mi estimate, post: regress `out' c.z_fm##i.geneticscore age physact i.ses i.sex z_ffm

esttab i_fm_`out' ///
       using `out'_fm_gene_int_pvalues_imputed.csv, replace ///
       cells("p(fmt(3))") ///
       label nodepvar ///
       varlabels(fm "FM") ///
	   keep(2.geneticscore#c.z_fm) ///
       nonumber alignment(D) noobs compress ///
	   plain
	   	   
}

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
	
eststo clear
eststo i_ffm_`out': quietly mi estimate, post: regress `out' c.z_ffm##i.geneticscore age physact i.ses i.sex z_fm

esttab i_ffm_`out' ///
       using `out'_ffm_gene_int_pvalues_imputed.csv, replace ///
       cells("p(fmt(3))") ///
       label nodepvar ///
       varlabels(ffm "FFM") ///
	   keep(2.geneticscore#c.z_ffm) ///
       nonumber alignment(D) noobs compress ///
	   plain
	   	   
}

//Create Combined Table//

preserve

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
import delimited "`out'_fm_gene_int_pvalues_imputed.csv", rowrange(3:3) colrange(2:) clear
rename i_fm_`out' p
save "`out'_fm_gene_int_pvalues_imputed.dta", replace
}

restore 

preserve

foreach out of varlist(lvm edv esv sv ef co hr lvws cmw csbp cpp dbp map svr hrvsdnn hrvrmsdd hrvti brs) {
import delimited "`out'_ffm_gene_int_pvalues_imputed.csv", rowrange(3:3) colrange(2:) clear
rename i_ffm_`out' p
save "`out'_ffm_gene_int_pvalues_imputed.dta", replace
}

restore

use "lvm_fm_gene_int_pvalues_imputed.dta", clear

append using "lvm_ffm_gene_int_pvalues_imputed.dta", force
append using "edv_fm_gene_int_pvalues_imputed.dta", force
append using "edv_ffm_gene_int_pvalues_imputed.dta", force
append using "esv_fm_gene_int_pvalues_imputed.dta", force
append using "esv_ffm_gene_int_pvalues_imputed.dta", force
append using "sv_fm_gene_int_pvalues_imputed.dta", force
append using "sv_ffm_gene_int_pvalues_imputed.dta", force
append using "ef_fm_gene_int_pvalues_imputed.dta", force
append using "ef_ffm_gene_int_pvalues_imputed.dta", force
append using "co_fm_gene_int_pvalues_imputed.dta", force
append using "co_ffm_gene_int_pvalues_imputed.dta", force
append using "hr_fm_gene_int_pvalues_imputed.dta", force
append using "hr_ffm_gene_int_pvalues_imputed.dta", force
append using "lvws_fm_gene_int_pvalues_imputed.dta", force
append using "lvws_ffm_gene_int_pvalues_imputed.dta", force
append using "cmw_fm_gene_int_pvalues_imputed.dta", force
append using "cmw_ffm_gene_int_pvalues_imputed.dta", force
append using "csbp_fm_gene_int_pvalues_imputed.dta", force
append using "csbp_ffm_gene_int_pvalues_imputed.dta", force
append using "cpp_fm_gene_int_pvalues_imputed.dta", force
append using "cpp_ffm_gene_int_pvalues_imputed.dta", force
append using "dbp_fm_gene_int_pvalues_imputed.dta", force
append using "dbp_ffm_gene_int_pvalues_imputed.dta", force
append using "map_fm_gene_int_pvalues_imputed.dta", force
append using "map_ffm_gene_int_pvalues_imputed.dta", force
append using "svr_fm_gene_int_pvalues_imputed.dta", force
append using "svr_ffm_gene_int_pvalues_imputed.dta", force
append using "hrvsdnn_fm_gene_int_pvalues_imputed.dta", force
append using "hrvsdnn_ffm_gene_int_pvalues_imputed.dta", force
append using "hrvrmsdd_fm_gene_int_pvalues_imputed.dta", force
append using "hrvrmsdd_ffm_gene_int_pvalues_imputed.dta", force
append using "hrvti_fm_gene_int_pvalues_imputed.dta", force
append using "hrvti_ffm_gene_int_pvalues_imputed.dta", force
append using "brs_fm_gene_int_pvalues_imputed.dta", force
append using "brs_ffm_gene_int_pvalues_imputed.dta", force

save "Supp_Table_4_int_pvalues.dta", replace
export delimited "Supp_Table_4_int_pvalues.csv", replace
