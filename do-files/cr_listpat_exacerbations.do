******************************************************
* EXACERBATIONS - CHRONIC DISEASES 
******************************************************

*Set gap: EVENTS WITHIN THIS GAP (DAYS) WILL COUNT AS ONE EXACERBATION
local gap = 14
global outcomes_chronic "bronch ild"
foreach outcome of global outcomes_chronic {

disp("Creating exacerbation episodes for outcome: `outcome'")

//Import data
use "$datafiles//cr_all_events_`outcome'_exacerbation.dta", clear 


// Remove duplicates 
duplicates drop

format %td obsdate 
gsort e_patid obsdate
desc
	
//Drop observations if record is from a hospitalisation as a secondary reason for admission 
cap drop if lrti_dx_type == 2 | pneumo_dx_type == 2

//Collapse a day per row for required vars
	if `"`outcome'"'  == "asthma" {
		cap drop if `outcome'_dx_type == 2 
		collapse (max) `outcome' `outcome'_dx_type `outcome'_exacerbation `outcome'_annual_review `outcome'_ocs lrti lrti_dx_type pneumo pneumo_dx_type ddate source, by(e_patid obsdate) 
	}

	if `"`outcome'"'  == "copd" {
		cap drop if `outcome'_dx_type == 2  
		collapse (max) copd copd_dx_type copd_exacerbation copd_annual_review copd_symptoms_breathlessness duration copd_symptoms_cough copd_symptoms_sputum copd_antibiotics copd_ocs lrti lrti_dx_type pneumo pneumo_dx_type ddate source, by(e_patid obsdate) 
	}	

	if `"`outcome'"'  == "bronch" {
		collapse (max) bronch bronch_dx_type bronch_exacerbation bronch_sec1  bronch_sec1_dx_type bronch_sec2 bronch_sec2_dx_type bronch_antibiotics duration lrti lrti_dx_type pneumo pneumo_dx_type ddate source, by(e_patid obsdate) 
	}

	if `"`outcome'"'  == "ild" {
		drop if `outcome'_dx_type ==2 
		collapse (max) ild ild_dx_type lrti lrti_dx_type pneumo pneumo_dx_type ddate source, by(e_patid obsdate) 	
	}

	list e_patid obsdate if e_patid == 3408816020460

//Generate end date (latest of discharge date (if hospitalised) and exacerbation gap)
gen gap = obsdate + `gap'
egen maxdate = rowmax(ddate gap)
format %td maxdate 

//Drop outcome specific annual review days
cap drop if `outcome'_annual_review == 1
cap drop `outcome'_annual_review 

	if `"`outcome'"'  == "asthma" {
		/*/ ALGORITHM:  ashtma exacerbation 
		if ON THE SAME DAY (not an annual review day), the patient recieved an: 
		1) `outcome' code (primary care, primary reason for hospitalisation or death) AND a prescription of oral corticosteroids, OR
		2) an exacerbation code in primary care, OR
		3) a code for lrti or pneumo (in primary care, primary reason for hospitalisation or death), OR
		4) was admitted to hospital with `outcome' as primary reason for hospitalisation, OR
		5) Death from `outcome'
		*/		
		keep if (`outcome' == 1 & `outcome'_ocs == 1) ///
			| (lrti ==1) | (pneumo ==1) ///
			| `outcome'_exacerbation == 1 /// 
			| (`outcome' == 1 & `outcome'_dx_type == 1) ///
			| (`outcome' == 1 & source == 4)
	}

	if `"`outcome'"'  == "copd" {
		/*/ ALGORITHM: copd exacerbation
		if ON THE SAME DAY (not an annual review day), the patient recieved: 
		1) a prescription for both antibiotics AND oral corticosteroids, OR
		2) 2 or more types of COPD symptoms in primary care and prescribed at least one copd treatment (antibiotics or oral corticosteroids), OR
		3) a copd exacerbation code, OR
		4) a code for lrti or pneumo (in primary care, primary reason for hospitalisation or death)
		5) was admitted to hospital with a hospitalisation with COPD as a primary diagnosis 
		6) Death from copd
		*/
		egen symptoms = rowtotal(copd_symptoms_breathlessness copd_symptoms_cough copd_symptoms_sputum)
		drop copd_symptoms_breathlessness copd_symptoms_cough copd_symptoms_sputum

		keep if (copd_antibiotics ==1 & copd_ocs == 1 & (duration >= 5 & duration <= 14)) ///
			| (symptoms >= 2 & (copd_antibiotics == 1 | copd_ocs == 1& (duration >= 5 & duration <= 14))) ///
			| copd_exacerbation == 1 ///
			| (lrti == 1) ///
			| (pneumo == 1) ///	
			| (copd == 1 & copd_dx_type == 1) ///
			| (copd == 1 & source == 4)
	}
		
	if `"`outcome'"'  == "bronch" {
		/*/ ALGORITHM: bronch exacerbation
		if ON THE SAME DAY, the patient recieved: 
		1) Bronchiectasis code AND antibiotic prescription; OR
		2) Bronchiectasis exacerbation code; OR
		3) a code for lrti or pneumo (in primary care, primary reason for hospitalisation or death)
		4) Admitted to hospital wih bronchiectasis in the primary position AND selected coded (sec2) codes in the secondary position
		5) Selected codes as the primary reason for diagnosis (sec1) and bronchiectasis in the second position
		6) Death from bronchiectasis
		*/
		keep if ///
			(bronch == 1 & (bronch_antibiotics == 1 & duration == 14)) ///
			| bronch_exacerbation == 1 ///
			| (lrti == 1) ///
			| (pneumo == 1) /// 
			| ((bronch == 1 & bronch_dx_type == 1) & (bronch_sec2 == 1 & bronch_sec2_dx_type == 2)) ///
			| ((bronch_sec1 == 1 & bronch_sec1_dx_type == 1) & (bronch == 1 & bronch_dx_type == 2)) ///
			| (bronch == 1 & source == 4)
	}	

	if `"`outcome'"' == "ild" {	
		/*/ ALGORITHM: ild exacerbation
		if ON THE SAME DAY, the patient recieved: 
		1) a code for lrti or pneumo (in primary care, primary reason for hospitalisation or death)
		2) ILD main reason for admission to hospital
		3) Death from ild
		*/
		keep if ///
			(lrti == 1) ///
			| (pneumo == 1) /// 
			| (ild == 1 & ild_dx_type == 1) ///
			| (ild == 1 & source == 4)
	}
	
//////////////////////////////////
*JOIN WITH COHORT PATIDS

desc 

tempfile temp
save `temp'

use "$cohort_patids_aurum", clear
joinby e_patid using `temp'

//Drop days where there is an event prior to indexdate (end of follow-up may be outcome specific, it is done later in R)
drop if obsdate < indexdate
////////////////////////////////////

		*prepare
		cap drop `outcome'_exc 
		cap drop `outcome'_exc_num 
		cap drop `outcome'_exc_contact 
		cap drop maxdate_lag
		cap drop date_diff_f

		gen `outcome'_exc = .
		gen `outcome'_exc_num = .
		gen `outcome'_exc_contact = .

		sort e_patid setid obsdate

		// maxdate of the previous row 
		bysort e_patid setid (obsdate): gen maxdate_lag = maxdate[_n-1]
		
		//obsdate of the previous row
	    bysort e_patid setid (obsdate): gen obsdate_lag = obsdate[_n-1]

		format %td maxdate_lag
		format %td obsdate_lag
				
	    list setid e_patid obsdate maxdate obsdate_lag maxdate_lag  if e_patid == 3408816020460  
	
		// drop rows where there is *FULL* overlap (obsdate to maxdate totally contained in a prior row)
		gen to_drop = 0 
		gen max_so_far = .
		
		*ensure dataset is sorted by obsdate
		sort e_patid setid (obsdate)
		
		*Generate running maximum
		bysort e_patid setid: replace max_so_far = max( max_so_far[_n-1], maxdate)
		
		*Drop rows that fully overlap
		replace to_drop = 1 if (obsdate > obsdate_lag) & (maxdate < max_so_far) 
		drop if to_drop == 1
		
		// recalculate maxdate 
		bysort e_patid setid (obsdate): replace maxdate_lag = maxdate[_n-1]
		
list setid e_patid obsdate maxdate obsdate_lag maxdate_lag max_so_far if e_patid == 3408816020460 
		
		//Group into one exacerbations if there is any overlap in time periods
	    bysort e_patid setid: replace `outcome'_exc = (maxdate_lag < obsdate) | missing(maxdate_lag)

		list setid e_patid obsdate maxdate obsdate_lag maxdate_lag `outcome'_exc  if e_patid == 3408816020460 
		
		// Calculate number (order) of exacerbations
		 bysort e_patid setid: replace `outcome'_exc_num = sum(`outcome'_exc)
		 
		 		list setid e_patid obsdate maxdate obsdate_lag maxdate_lag `outcome'_exc `outcome'_exc_num if e_patid == 3408816020460 
		 
       // Calculate number of contacts within an exacerbation
		bysort e_patid setid `outcome'_exc_num: replace `outcome'_exc_contact = _N

		disp("Number of rows to be grouped")

		list setid e_patid obsdate maxdate `outcome'_exc `outcome'_exc_num `outcome'_exc_contact if e_patid == 3408816020460 
		
		// Collapse to one row per exacerbation
		collapse (min) obsdate (max) `outcome'_exc maxdate `outcome'_exc_contact, by(e_patid setid indexdate exposed `outcome'_exc_num) 
		
		list setid e_patid obsdate maxdate `outcome'_exc `outcome'_exc_num `outcome'_exc_contact if e_patid == 3408816020460 

//	}

//See resulting dataset
desc
tab `outcome'_exc
	
save "$datafiles//cr_listpat_`outcome'_exacerbations.dta", replace

bysort e_patid setid (`outcome'_exc_num): gen maxdate_lag = maxdate[_n-1]

gen date_diff = obsdate - maxdate_lag
list if date_diff <=0 

count if date_diff <0
count if date_diff <= 0 

drop date_diff
drop maxdate_lag

}


clear