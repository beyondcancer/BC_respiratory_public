/********************************************************************
APPEND ALL AURUM, RX, HES & ONS EVENTS INTO ONE DATASET FOR EACH RESPIRATORY
OUTCOME
********************************************************************/
foreach outcome of global all_codelists {
	

    di "Creating cr_all events file for outcome: `outcome'"

    local dataset_found = 0

    // Attempt to use the dx dataset
    capture confirm file "$datafiles\cr_all_events_dx_`outcome'_aurum.dta"
    if _rc == 0 {
        di "File exists: cr_all_events_dx_`outcome'_aurum.dta"
        use "$datafiles\cr_all_events_dx_`outcome'_aurum.dta", clear
        gen source = 1
        destring e_patid, replace
        local dataset_found = 1
        local initial_dataset = "dx"
    } 
    else {
        di "File not found: cr_all_events_dx_`outcome'_aurum.dta"
    }

    // If dx doesn't exist, attempt to use the rx dataset
    if `dataset_found' == 0 {
        capture confirm file "$datafiles\cr_all_events_rx_`outcome'_aurum.dta"
        if _rc == 0 {
            di "File exists: cr_all_events_rx_`outcome'_aurum.dta"
            use "$datafiles\cr_all_events_rx_`outcome'_aurum.dta", clear
            gen source = 2
            local dataset_found = 1
            local initial_dataset = "rx"
        } 
        else {
            di "File not found: cr_all_events_rx_`outcome'_aurum.dta"
        }
    }

    // If rx doesn't exist, attempt to use the hes dataset
    if `dataset_found' == 0 {
        capture confirm file "$datafiles\cr_all_events_hes_`outcome'_aurum.dta"
        if _rc == 0 {
            di "File exists: cr_all_events_hes_`outcome'_aurum.dta"
            use "$datafiles\cr_all_events_hes_`outcome'_aurum.dta", clear
            gen source = 3
            local dataset_found = 1
            local initial_dataset = "hes"
        } 
        else {
            di "File not found: cr_all_events_hes_`outcome'_aurum.dta"
        }
    }

    // If hes doesn't exist, attempt to use the ons dataset
    if `dataset_found' == 0 {
        capture confirm file "$datafiles\cr_all_events_ons_`outcome'_aurum.dta"
        if _rc == 0 {
            di "File exists: cr_all_events_ons_`outcome'_aurum.dta"
            use "$datafiles\cr_all_events_ons_`outcome'_aurum.dta", clear
            gen source = 4
            local dataset_found = 1
            local initial_dataset = "ons"
        } 
        else {
            di "File not found: cr_all_events_ons_`outcome'_aurum.dta"
        }
    }

    // If none of the files exist, skip to the next iteration
    if `dataset_found' == 0 {
        di "Warning: No available dataset found for outcome: `outcome'"
        continue
    }

    // List of all datasets to potentially append
    local datasets_to_append = "rx hes ons"

    foreach dataset in `datasets_to_append' {
		
        // Skip appending the initial dataset
        if "`dataset'" != "`initial_dataset'" {
            capture confirm file "$datafiles\cr_all_events_`dataset'_`outcome'_aurum.dta"
            if _rc == 0 {
                di "Appending file: cr_all_events_`dataset'_`outcome'_aurum.dta"
                append using "$datafiles\cr_all_events_`dataset'_`outcome'_aurum.dta"
                // Set appropriate source value for each dataset
                if "`dataset'" == "rx" {
                    recode source .= 2
                }
                else if "`dataset'" == "hes" {
                    recode source .= 3
                }
                else if "`dataset'" == "ons" {
                    recode source .= 4
                }
            } 
            else {
                di "File not found, skipping: cr_all_events_`dataset'_`outcome'_aurum.dta"
            }
        }
    }

    // Define labels for the source variable if they exist
    cap lab define source 1 dx_aurum 2 rx 3 hes 4 ons
    lab val source source

// Check if ddate and duration exist

capture confirm variable duration
local has_duration = _rc == 0

capture confirm variable out_cat
local has_outcat = _rc == 0

capture confirm variable ddate
local has_ddate = _rc == 0

// Construct keep statement based on existence of variables
local keepvars e_patid obsdate `outcome'* source

if `has_duration' local keepvars `keepvars' duration 
if `has_outcat' local keepvars `keepvars'  out_cat certlevel
if `has_ddate' local keepvars `keepvars' ddate 

keep `keepvars'
order `keepvars'

    // Drop missing data and duplicates
    cap destring e_patid, replace
	drop if e_patid == .
    drop if obsdate == .
    duplicates drop

    // Describe and list data
    desc
    list in 1/10

    // Save the final dataset
    save "$datafiles/cr_all_events_`outcome'_aurum.dta", replace
}



/********************************************************************
APPEND ALL AURUM, RX, HES & ONS EVENTS INTO ONE DATASET FOR EACH RESPIRATORY
OUTCOME EXACERBATION
********************************************************************/


foreach outcome of global outcomes_chronic {
    use "$datafiles//cr_all_events_`outcome'_aurum.dta", clear

    // Append rows from other datasets
	*annual review
    capture confirm file "$datafiles//cr_all_events_`outcome'_annual_review_aurum.dta"
    if _rc == 0 {
        append using  "$datafiles//cr_all_events_`outcome'_annual_review_aurum.dta"
		desc
    } 
	else {
        di "File not found: cr_all_events_`outcome'_annual_review_aurum.dta"
    }
	*exacerbation
    capture confirm file "$datafiles//cr_all_events_`outcome'_exacerbation_aurum.dta"
    if _rc == 0 {
        append using  "$datafiles//cr_all_events_`outcome'_exacerbation_aurum.dta"
		desc
	else {
        di "File not found: cr_all_events_`outcome'_exacerbation_aurum.dta"
    }
	*symptoms
	    capture confirm file "$datafiles//cr_all_events_`outcome'_symptoms_aurum.dta"
    if _rc == 0 {
        append using     "$datafiles//cr_all_events_`outcome'_symptoms_aurum.dta"
		desc
    } 
	else {
        di "File not found: cr_all_events_`outcome'_symptoms_aurum.dta"
    }
	*antibiotics
	    capture confirm file "$datafiles//cr_all_events_`outcome'_antibiotics_aurum.dta"
    if _rc == 0 {
        append using "$datafiles//cr_all_events_`outcome'_antibiotics_aurum.dta"
		desc
    } 
	else {
        di "File not found: cr_all_events_`outcome'_antibiotics_aurum.dta"
    }
	*oral corticosteroids
	    capture confirm file "$datafiles//cr_all_events_`outcome'_ocs_aurum.dta"
    if _rc == 0 {
        append using "$datafiles//cr_all_events_`outcome'_ocs_aurum.dta"
		desc
    } 
	else {
        di "File not found: cr_all_events_`outcome'_ocs_aurum.dta"
    }
	*secondary codes (bronch only)
	    capture confirm file "$datafiles//cr_all_events_hes_`outcome'_sec2_aurum.dta"
    if _rc == 0 {
        append using "$datafiles//cr_all_events_hes_`outcome'_sec2_aurum.dta"
		desc
    } 
	else {
        di "File not found: cr_all_events_hes_`outcome'_sec2_aurum"
    }
	    capture confirm file "$datafiles//cr_all_events_hes_`outcome'_sec1_aurum.dta"
    if _rc == 0 {
        append using "$datafiles//cr_all_events_hes_`outcome'_sec1_aurum.dta"
		desc
    } 
	else {
        di "File not found: cr_all_events_hes_`outcome'_sec1_aurum.dta"
    }
	*lower respiratory tract infections
    capture confirm file "$datafiles//cr_all_events_lrti_aurum.dta"
    if _rc == 0 {
        append using "$datafiles//cr_all_events_lrti_aurum.dta"
		desc
    } 
	else {
        di "File not found: cr_all_events_lrti_aurum.dta"
    }
	*pneumonia
    capture confirm file "$datafiles//cr_all_events_pneumo_aurum.dta"
    if _rc == 0 {
        append using "$datafiles//cr_all_events_pneumo_aurum.dta"
		desc
    } 
	else {
        di "File not found: cr_all_events_pneumo_aurum.dta"
    }

//save the final dataset
save "$datafiles//cr_all_events_`outcome'_exacerbation.dta", replace

}
}

clear
/*
/********************************************************************
ERASE INTERMEDIATE FILES
********************************************************************/

foreach outcome of global aurum_codelists {
local dx_out : dir "$datafiles" files "cr_all_events_dx_`outcome'_aurum*.dta"
foreach file of local dx_out {
    display "Removing file: `file'"   
    erase "$datafiles//`file'"
}
}

foreach outcome of global rx_codelists {
local dx_out : dir "$datafiles" files "cr_all_events_rx_`outcome'_aurum*.dta"
foreach file of local dx_out {
    display "Removing file: `file'"   
    erase "$datafiles//`file'"
}
}

foreach outcome of global hes_codelists {
local hes_out : dir "$datafiles" files "cr_all_events_hes*_`outcome'_aurum*.dta"
foreach file of local hes_out {
    display "Removing file: `file'"   
    erase "$datafiles//`file'"
}
}

foreach outcome of global ons_codelists {
local ons_out : dir "$datafiles" files"cr_all_events_ons_`outcome'*.dta"
foreach file of local ons_out {
    display "Removing file: `file'"   
    erase "$datafiles//`file'"
}
}


log close
	