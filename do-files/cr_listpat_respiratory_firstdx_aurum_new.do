clear
capture log close
log using "$logfiles_create_dataset\cr_listpat_respiratory_aurum.txt", replace text

/************************************************************
Outcomes FIRST_EVENT analysis:

DIAGNOSIS:  
	- outcome_inc_dx: at least one relevant code denoting incident respiratory condition after cancer in CPRD :At least one relevant code in any diagnosis field in CPRD Aurum or HES-APC. 

	- outcome_hosp: At least one code denoting a hospitalisation with the respiratory condition after cancer as the primary reason for admittance that occurs after index date
	 
	- outcome_death: At least one code with the outcome as a cause of death
	
*EXCLUSION flags: 
- Flag those with "history of" (b_) or with a prevalence code (denoted by variable out_cat: dx_prev) as the first code after index date (var: first_event_outcome_type)
	
********************************************************************/


foreach outcome of global outcomes_all {
	
di "Create listpat for outcome: `outcome'"

frame change default

use "$cohort_patids_aurum", clear

keep e_patid setid indexdate exposed
count
joinby e_patid using "$datafiles\cr_all_events_`outcome'_aurum"

tab exposed `outcome'

format obsdate %td "DDmmmYYY"

********************************************************************
*DIAGNOSIS - includes primary care, hosp and ons codes - ONLY 1st diagnosis 
********************************************************************
di "Creating listpat with first dx (one row per patient): `outcome'"

// Generate flag for events that occur before the indexdate
if strpos("$outcomes_chronic", "`outcome'") > 0 {  
gen b_`outcome'= 1 if obsdate < indexdate & `outcome'_dx_type != 2
}

if strpos("$outcomes_infection", "`outcome'") > 0 {  
gen b_`outcome'= 1 if obsdate < indexdate & source != 1
}

gen dol_b_`outcome' = obsdate if b_`outcome' == 1 
format dol_b_`outcome' %td "DDmmmYYY"

// Generate flag for events that occur after index date and do not denote prevalent disease (can be primary care, hospitalisation or death).
if strpos("$outcomes_chronic", "`outcome'") > 0 {  
gen any_`outcome'_inc_dx = 1 if obsdate >= indexdate & out_cat != 1 & `outcome'_dx_type != 2
}
if strpos("$outcomes_infection", "`outcome'") > 0 { 
gen any_`outcome'_inc_dx = 1 if obsdate >= indexdate & out_cat != 1 & source != 1
}

gen dof_any_`outcome'_inc_dx = obsdate if any_`outcome'_inc_dx==1 
format dof_any_`outcome'_inc_dx  %td "DDmmmYYY"

// Generate flag for events that occur after index date in primary care and do not denote prevalent disease (can be primary care, hospitalisation or death) - NOT FOR INFECTIONS

if strpos("$outcomes_chronic", "`outcome'") > 0 {  
gen `outcome'_inc_dx = 1 if obsdate >= indexdate & out_cat != 1 & source == 1
gen dof_`outcome'_inc_dx = obsdate if `outcome'_inc_dx==1 
format dof_`outcome'_inc_dx  %td "DDmmmYYY"
}

if strpos("$outcomes_infection", "`outcome'") > 0 {  
gen `outcome'_inc_dx = . 
gen dof_`outcome'_inc_dx = .
format dof_`outcome'_inc_dx  %td "DDmmmYYY"
}

// Generate flag for hospitalisation events that occur after index date and do not denote prevalent disease.
if strpos("$outcomes_chronic", "`outcome'") > 0 {  
gen `outcome'_hosp = 1 if obsdate >= indexdate & source == 3 & `outcome'_dx_type == 1
}

if strpos("$outcomes_infection", "`outcome'") > 0 { 
gen `outcome'_hosp = 1 if obsdate >= indexdate & source == 3
}

gen dof_`outcome'_hosp = obsdate if `outcome'_hosp==1 
format dof_`outcome'_hosp  %td "DDmmmYYY"


// Generate flag for death events that occur after index date and do not denote prevalent disease.
gen `outcome'_death = 1 if obsdate >= indexdate & source == 4 
gen dof_`outcome'_death = obsdate if `outcome'_death==1 
format dof_`outcome'_death  %td "DDmmmYYY"

// Generate flag for prevalent dx_codes
gen `outcome'_prev = 1 if obsdate >= indexdate & out_cat == 1 
gen dof_`outcome'_prev = obsdate if `outcome'_prev==1
format dof_`outcome'_prev  %td "DDmmmYYY"

// Generate the type for the first ever dx 
bysort e_patid setid (obsdate): gen first_`outcome'_cat = out_cat if _n==1 

// Identify the first event (after index date) type based on the minimum date
egen first_event_date = rowmin(dof_`outcome'_inc_dx dof_`outcome'_hosp dof_`outcome'_prev dof_`outcome'_death)
gen first_event_`outcome'_type = ""
replace first_event_`outcome'_type = "inc_dx" if first_event_date == dof_`outcome'_inc_dx & b_`outcome' != 1
replace first_event_`outcome'_type = "hosp" if first_event_date == dof_`outcome'_hosp & b_`outcome' != 1
replace first_event_`outcome'_type = "death" if first_event_date == dof_`outcome'_death & b_`outcome' != 1
replace first_event_`outcome'_type = "prev" if first_event_date == dof_`outcome'_prev & b_`outcome' != 1

format first_event_date %td "DDmmmYYY"

sort e_patid setid (obsdate)

//select firstdx after index date
collapse (max) b_`outcome' any_`outcome'_inc_dx `outcome'_prev `outcome'_inc_dx `outcome'_hosp `outcome'_death dol_b_`outcome' (min) dof_any_`outcome'_inc_dx dof_`outcome'_prev dof_`outcome'_inc_dx dof_`outcome'_hosp dof_`outcome'_death first_`outcome'_cat (firstnm) first_event_`outcome'_type, by(e_patid setid indexdate exposed)			

list in 1/5

drop if dof_any_`outcome'_inc_dx > date("29mar2021", "DMY")

hist dof_any_`outcome'_inc_dx

order setid e_patid exposed indexdate b_`outcome' dol_b_`outcome' `outcome'_prev `outcome'_inc_dx dof_`outcome'_inc_dx `outcome'_hosp dof_`outcome'_hosp `outcome'_death dof_`outcome'_death first_`outcome'_cat first_event_`outcome'_type

keep setid e_patid exposed indexdate b_`outcome' dol_b_`outcome' any_`outcome'_inc_dx dof_any_`outcome'_inc_dx first_event_`outcome'_type
order setid e_patid exposed indexdate b_`outcome' dol_b_`outcome' any_`outcome'_inc_dx dof_any_`outcome'_inc_dx first_event_`outcome'_type

desc
tab any_`outcome'_inc_dx first_event_`outcome'_type, m

save "$datafiles/cr_listpat_`outcome'_aurum_firstdx", replace
}

log close 
clear
