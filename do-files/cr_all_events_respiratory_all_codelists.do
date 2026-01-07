clear
capture log close
log using "$logfiles_create_dataset\cr_all_events_respiratory_aurum.txt", replace text


/********************************************************************
CREATE FILE WITH ALL RESPIRATORY OUTCOME EVENTS (inc definite, probable, possible, h/o)
 IN aurum, HES and ONS
********************************************************************/

/********************************************************************
CPRD AURUM Dx for outcomes - note: Only for chronic outcomes
********************************************************************/

*1. Append all Dx covariate codelists creating variable with codelist name: aurum

/*Convert codelist to lookup with numeric codes*/

	foreach outcome of global aurum_codelists {
		di "creating codelist with lshtmcode for `outcome'"
		use "$codelists//cr_codelist_`outcome'_aurum.dta", clear 				/* HC = where do these codelists come from? */
		tostring medcodeid, force replace
		cap drop `outcome'
		gen `outcome' = 1
		//rename cat_resp `outcome'_cat
		duplicates drop medcodeid, force
		merge 1:1 medcodeid using "$datafiles//aurummedicaldict_lookup.dta", keep(match) nogen
		list in 1
		save "$codelists//cr_codelist_`outcome'_aurum_lshtmcode", replace
		clear	
	} // outcomes

/*Remove files from directory before starting*/

	foreach outcome of global aurum_codelists {
	local dx_out : dir "$datafiles" files "cr_all_events_dx_`outcome'_aurum*.dta"
	foreach file of local dx_out {
		display "Removing file: `file'"   
		cap erase "`file'"
	}
	}

/* Search all obs files (rawdata) */

	foreach outcome of global aurum_codelists  {
	di "searching primary care files for: `outcome' dx events"
		foreach part in 1 2 3 4 5 6 7 8 9 10 {
				foreach file of numlist 1/267 {
					di "Part `part ', Observation file `file'"
					quietly {
						if `file'<10 {
							use "$rawdata_aurum\\part_`part'\\Observation\\part`part'_Observation_00`file'.dta"
							}
							
						if `file'>9 & `file'<100 {
							use "$rawdata_aurum\\part_`part'\\Observation\\part`part'_Observation_0`file'.dta"
							}
							
						if `file'>99 {
							cap use "$rawdata_aurum\\part_`part'\\Observation\\part`part'_Observation_`file'.dta"
							if _rc==601 {
							continue
							}				
							}
						merge m:1 lshtmcode using  "$codelists//cr_codelist_`outcome'_aurum_lshtmcode", keep(match) nogen
						save "$datafiles//cr_all_events_dx_`outcome'_aurum_`part'_`file'.dta", replace		
						
					} // quietly
				} // part
			} // file
				
		********************************************************************************
		* Create a single file with all events
		********************************************************************************
		/* output files */
		clear
		set obs 1
		gen e_patid=""
		gen obsdate=.
		gen `outcome'=.
		save "$datafiles//cr_all_events_dx_`outcome'_aurum_1.dta", replace
		save "$datafiles//cr_all_events_dx_`outcome'_aurum_2.dta", replace
		save "$datafiles//cr_all_events_dx_`outcome'_aurum_3.dta", replace
		save "$datafiles//cr_all_events_dx_`outcome'_aurum_4.dta", replace
		save "$datafiles//cr_all_events_dx_`outcome'_aurum_5.dta", replace
		save "$datafiles//cr_all_events_dx_`outcome'_aurum_6.dta", replace
		save "$datafiles//cr_all_events_dx_`outcome'_aurum_7.dta", replace
		save "$datafiles//cr_all_events_dx_`outcome'_aurum_8.dta", replace
		save "$datafiles//cr_all_events_dx_`outcome'_aurum_9.dta", replace
		save "$datafiles//cr_all_events_dx_`outcome'_aurum_10.dta", replace

		/* append the different files in each part */
		foreach part in 1 2 3 4 5 6 7 8 9 10 {
			use "$datafiles//cr_all_events_dx_`outcome'_aurum_`part'.dta", clear
					
				foreach file of numlist 1/267 {
				cap append using  "$datafiles//cr_all_events_dx_`outcome'_aurum_`part'_`file'.dta"
					if _rc==601 {
					continue
					}
				}
			save "$datafiles//cr_all_events_dx_`outcome'_aurum_`part'.dta", replace
			}	
				
		/* append all parts together */	
		use "$datafiles//cr_all_events_dx_`outcome'_aurum_1.dta", clear
		
			foreach part in 2 3 4 5 6 7 8 9 10 {
			append using "$datafiles//cr_all_events_dx_`outcome'_aurum_`part'.dta"
			}


		/* SAVE FILE WITH ALL EVENTS */
		duplicates drop
		drop if obsdate==.
		drop if e_patid==""
		drop lshtmcode
		keep e_patid obsdate `outcome'* out_cat certlevel 
		order e_patid obsdate `outcome'* out_cat certlevel 
		compress
		save "$datafiles//cr_all_events_dx_`outcome'_aurum.dta", replace


		/* Erase intermediate files*/
			foreach part in 1 2 3 4 5 6 7 8 9 10 {
				foreach file of numlist 1/267 {
					cap erase "$datafiles//cr_all_events_dx_`outcome'_aurum_`part'_`file'.dta"
				}	
			}
			
			foreach part in 1 2 3 4 5 6 7 8 9 10 {
				cap erase "$datafiles//cr_all_events_dx_`outcome'_aurum_`part'.dta"
			}	
		
	} //outcome


/********************************************************************
HES APC Dx for ALL RESPIRATORY OUTCOMES

Interested in:
the primary discharge Dx (from primary diag hosp) 
and other Dx in hospital episode (from hes_epi)
********************************************************************/

/*generate outcome var for each codelists*/
	foreach outcome of global hes_codelists  {
	di "preparing hes codelists for: `outcome'"
	use "$codelists\cr_codelist_`outcome'_hes.dta", clear
	cap gen `outcome' = 1
	cap rename cat_resp `outcome'_cat
	duplicates drop
	save "$codelists\cr_codelist_`outcome'_hes.dta", replace
	}
	
/*remove files before starting*/
	foreach outcome of global hes_codelists  {
	local hes_out : dir "$datafiles" files "cr_all_events_hes*_`outcome'_aurum*.dta"
		foreach file of local hes_out {
			display "Removing file: `file'"   
			cap erase "`file'"
		}
	}
	
********************************************************************************
* Look for events in primary dx file
*******************************************************************************

di "HES events"

/* load in primary dx file */		
use "$rawdata_aurum_linked\e_aurum_hes_primary_diag_hosp_20_000268_dm.dta", clear
rename icd_primary icd
rename admidate obsdate

sort e_patid (obsdate) 

gen ddate= date(discharged, "DMY")
format ddate %td "DDMMMYYYY"

list if obsdate != ddate in 1/10

/* create new frame to execute the merge */		
cap frame drop merge_frame

frame create merge_frame

/*loop through all outcomes*/	

foreach outcome of global hes_codelists  {
	
	/*Make sure the loop is using the merge frame*/		
	frame change merge_frame

	/*load in the hes epi data into the merge_frame*/	
	frame copy default merge_frame, replace
	
	di "checking for diagnosis as a primary diagnosis in HES for `outcome':"	
	
	/*perform merge*/
	merge m:1 icd using "$codelists\cr_codelist_`outcome'_hes.dta", keep(match) nogen
	
	/*create indicator variable for outcome that is from primary diagnosis file*/
	
	gen `outcome'_dx_type = 1 // indicates outcome is from primary diagnosis file
	
	/*keep relevant vars and save file*/
	keep e_patid obsdate ddate `outcome'* out_cat certlevel icd 
	drop if obsdate==.
	drop if e_patid==.
	save "$datafiles//cr_all_events_hes_pri_`outcome'_aurum.dta", replace
	desc
	
	/*change back to the main frame*/
	frame change default
} //outcome

********************************************************************************
* Look for events in epi dx file
*******************************************************************************

/*make sure the working frame is the default frame*/
assert "`c(frame)'" == "default"

di "HES EPI"
/*loop through all outcomes*/
foreach outcome of global hes_sec_codelists {
	

/* load in epi dx file */

	use "$rawdata_aurum_linked\e_aurum_hes_diagnosis_epi_20_000268_dm.dta", clear
	rename epistart obsdate
	keep e_patid icd d_order obsdate epiend

	sort e_patid (obsdate) epiend	
														 /* HC added: you need to sort the data by event date as well */ 
	bysort e_patid (obsdate): gen ddate = epiend[_N]  
	*format ddate %tddd "DDMMMYYYY"
	format ddate %td

	di "checking for diagnosis in any discharge field in HES for `outcome':"	
	
	/*perform merge*/
	merge m:1 icd using "$codelists\cr_codelist_`outcome'_hes.dta", keep(match) nogen
	
	/*create indicator variable for outcome that is from primary diagnosis file*/
	gen `outcome'_dx_type = 2 
	desc

	/*keep relevant vars and save file*/
	keep e_patid obsdate ddate `outcome'* out_cat certlevel icd d_order 
	drop if obsdate==.
	drop if e_patid==.
	drop if d_order == 1
	save "$datafiles//cr_all_events_hes_epi_`outcome'_aurum.dta", replace

	clear 
}

/*make sure the working frame is the default frame*/
assert "`c(frame)'" == "default"


********************************************************************************
* APPEND HES FILES TOGETHER
********************************************************************************


foreach outcome of global hes_codelists {										/* HC = 1) not sure why you would need this code, it is impossible to not have events for the outcome in HES with such a big sample. 2) Lines 297 and 303 are not independent conditions of line 286. This means that you file will be replaced by the latest condition that appears. This needs to be checked. */
    
    di "Processing HES datasets for outcome: `outcome'"
    
    // Check if the epi file exists
    capture confirm file "$datafiles//cr_all_events_hes_epi_`outcome'_aurum.dta"
    local epi_exists = _rc == 0
    
    // Check if the pri file exists
    capture confirm file "$datafiles//cr_all_events_hes_pri_`outcome'_aurum.dta"
    local pri_exists = _rc == 0
    
    if `epi_exists' & `pri_exists' {
        di "Both files exist for outcome: `outcome'. Appending epi to pri."
        
        use "$datafiles//cr_all_events_hes_pri_`outcome'_aurum.dta", clear
        append using "$datafiles//cr_all_events_hes_epi_`outcome'_aurum.dta"
        
        keep e_patid obsdate ddate `outcome'* out_cat certlevel 
        order e_patid obsdate ddate `outcome'* out_cat certlevel 
        
        save "$datafiles//cr_all_events_hes_`outcome'_aurum.dta", replace
    } 
	
    else if `pri_exists' & !`epi_exists' {
        di "Only pri file exists for outcome: `outcome'. Saving pri file as new dataset."
        
        use "$datafiles//cr_all_events_hes_pri_`outcome'_aurum.dta", clear
        save "$datafiles//cr_all_events_hes_`outcome'_aurum.dta", replace
    } 
    else if `epi_exists' & !`pri_exists' {
        di "Only epi file exists for outcome: `outcome'. Saving epi file as new dataset."
        
        use "$datafiles//cr_all_events_hes_epi_`outcome'_aurum.dta", clear
        save "$datafiles//cr_all_events_hes_`outcome'_aurum.dta", replace
    } 
    else {
        di "Neither file exists for outcome: `outcome'. Skipping processing."
    }
}


********************************************************************
****CPRD Rx
********************************************************************
* 3. Append all Rx covariate codelists creating variable with codelist name

foreach outcome of global rx_codelists  {
*  create temp file
    clear
	set obs 0
	gen prodcodeid = .
    save "$codelists/temp_prodcodelist_aurum.dta", replace

    use "$codelists/cr_codelist_`outcome'_rx.dta", clear
	desc
	gen `outcome' = 1
	keep prodcodeid `outcome'

    * Convert prodcodeid to numeric if it's not already
    cap destring prodcodeid, replace

    * Append to the temporary file
    append using "$codelists/temp_prodcodelist_aurum.dta"

    * Save the updated temporary file
    save "$codelists/temp_prodcodelist_aurum.dta", replace

distinct prodcodeid

save "$datafiles/cr_rx_`outcome'_aurum.dta", replace  

* convert codelist to lookup with numeric codes
use "$datafiles\cr_rx_`outcome'_aurum.dta", clear
tostring prodcodeid, replace format("%20.0f")
merge m:1 prodcodeid using "$datafiles//aurumproddict_lookup.dta", keep(match) nogen
destring prodcodeid, replace
save "$codelists\cr_rx_`outcome'_aurum_lshtmcode.dta", replace

//this makes the final list of prescription codes for the respiratory conditions of interest only

********************************************************************************
clear
* search all obs files
foreach part in 1 2 3 4 5 6 7 8 9 10 {
	foreach file of numlist 1/249 {
	
		di "Part `part ', DrugIssue file `file'"
		
		quietly {
			if `file'<10 {
				use "$rawdata_aurum\\part_`part'\\DrugIssue\\part`part'_DrugIssue_00`file'.dta"
				}
				
			if `file'>9 & `file'<100 {
				use "$rawdata_aurum\\part_`part'\\DrugIssue\\part`part'_DrugIssue_0`file'.dta"
				}
				
			if `file'>99 {
				cap use "$rawdata_aurum\\part_`part'\\DrugIssue\\part`part'_DrugIssue_`file'.dta"
				if _rc==601 {
				continue
				}				
				}
				
			joinby lshtpcode using  "$codelists//cr_rx_`outcome'_aurum_lshtmcode"
			sum
			save "$datafiles//cr_rx_`outcome'_aurum_`part'_`file'.dta", replace			
			} // quietly			
	} // part
	} // file
	
********************************************************************************
* output file
clear
set obs 1
gen e_patid=""
gen obsdate=.
gen lshtpcode=.
gen `outcome' = .
save "$datafiles//cr_rx_`outcome'_aurum_1.dta", replace
save "$datafiles//cr_rx_`outcome'_aurum_2.dta", replace
save "$datafiles//cr_rx_`outcome'_aurum_3.dta", replace
save "$datafiles//cr_rx_`outcome'_aurum_4.dta", replace
save "$datafiles//cr_rx_`outcome'_aurum_5.dta", replace
save "$datafiles//cr_rx_`outcome'_aurum_6.dta", replace
save "$datafiles//cr_rx_`outcome'_aurum_7.dta", replace
save "$datafiles//cr_rx_`outcome'_aurum_8.dta", replace
save "$datafiles//cr_rx_`outcome'_aurum_9.dta", replace
save "$datafiles//cr_rx_`outcome'_aurum_10.dta", replace


foreach part in 1 2 3 4 5 6 7 8 9 10 {
	use "$datafiles//cr_rx_`outcome'_aurum_`part'.dta"

	foreach file of numlist 1/249 {
		cap append using  "$datafiles//cr_rx_`outcome'_aurum_`part'_`file'.dta"
		if _rc==601 {
			continue
			}
		}
	save "$datafiles//cr_rx_`outcome'_aurum_`part'.dta", replace
}	
	
	
use "$datafiles//cr_rx_`outcome'_aurum_1.dta"
forvalues var = 2/10 {
append using "$datafiles//cr_rx_`outcome'_aurum_`var'.dta"
}

destring e_patid, replace
keep e_patid obsdate duration lshtpcode `outcome'
count

save "$datafiles\cr_all_events_rx_`outcome'_aurum.dta", replace

********************************************************************************
* erase intermediate files
foreach part in 1 2 3 4 5 6 7 8 9 10 {
	foreach file of numlist 1/267 {
		cap erase "$datafiles//cr_rx_`outcome'_aurum_`part'_`file'.dta"
}	
}

foreach part in 1 2 3 4 5 6 7 8 9 10 {
	cap erase "$datafiles//cr_rx_`outcome'_aurum_`part'.dta"
}
	
}


/********************************************************************
ONS DX for ALL RESPIRATORY OUTCOMES
********************************************************************/
foreach outcome of global ons_codelists  {
	local ons_out : dir "$datafiles" files "cr_all_events_ons_`outcome'*.dta"
	foreach file of local ons_out {
		display "Removing file: `file'"   
		cap erase "`file'"
	}
}

/*loop through all outcomes*/	
foreach outcome of global ons_codelists   {
	
	/*generate outcome var for each codelists*/
	di "preparing ONS codelists for outcome `outcome' "
	use "$codelists\cr_codelist_`outcome'_ons.dta", clear
	gen `outcome' = 1
	rename cat_resp `outcome'_cat
	duplicates drop icd, force 
	
	di "looking for all deaths with `outcome' as a cause"
	
	/*Merge with ONS data*/
	rename icd cause
	joinby cause using "$rawdata_aurum_linked\e_aurum_death_patient_20_000268_dm.dta"
	rename dod obsdate
	rename cause icd

	/*keep relevant vars and save file*/
	keep e_patid obsdate `outcome'* out_cat certlevel 
	order e_patid obsdate `outcome'* out_cat certlevel 
	drop if obsdate==.
	drop if e_patid==.

	desc
	save "$datafiles\cr_all_events_ons_`outcome'_aurum.dta", replace
	} // outcome

/*	
/********************************************************************
APPEND ALL AURUM, RX, HES & ONS EVENTS INTO ONE DATASET FOR EACH RESPIRATORY
OUTCOME
********************************************************************/
foreach outcome of global all_codelists {

di "creating cr_all events file for outcome : `outcome' "


    capture confirm file "$datafiles\cr_all_events_dx_`outcome'_aurum.dta"

    if _rc == 0 {
	
use "$datafiles\cr_all_events_dx_`outcome'_aurum.dta", clear
gen source =1
destring e_patid, replace

append using "$datafiles\cr_all_events_hes_`outcome'_aurum.dta"
recode source .=2
	} 
	
	else {
		
		capture confirm file "$datafiles\cr_all_events_rx_`outcome'_aurum.dta"

    if _rc == 0 {
		
use "$datafiles\cr_all_events_hes_`outcome'_aurum.dta", clear
gen source=2

	}
	
	else {
		
append using "$datafiles\cr_all_events_ons_`outcome'_aurum.dta"
recode source .=3

cap lab define source 1 dx_aurum 2 hes 3 ons
lab val source source

keep e_patid obsdate ddate `outcome'*  out_cat source certlevel
order e_patid obsdate ddate `outcome'* out_cat source certlevel

drop if e_patid == .
drop if obsdate == .
duplicates drop

desc
list in 1/10
save "$datafiles//cr_all_events_`outcome'_aurum.dta", replace
}

/********************************************************************
ERASE INTERMEDIATE FILES
********************************************************************/

foreach outcome of local outcomes_chronic {
local dx_out : dir "$datafiles" files "cr_all_events_dx_`outcome'_aurum*.dta"
foreach file of local dx_out {
    display "Removing file: `file'"   
   // erase "`file'"
}
}

foreach outcome of local outcomes_chronic {
local hes_out : dir "$datafiles" files "cr_all_events_hes*_`outcome'_aurum*.dta"
foreach file of local hes_out {
    display "Removing file: `file'"   
   // erase "`file'"
}

local ons_out : dir "$datafiles" files"cr_all_events_ons_`outcome'*.dta"
foreach file of local ons_out {
    display "Removing file: `file'"   
   // erase "`file'"
}
}


log close
	
	