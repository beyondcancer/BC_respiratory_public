
********************************************************************************
* do file author:	Kirsty Andresen
* Date: 			19 September 2024
* Description: 		Master file for BC_respiratory_outcomes in stata
********************************************************************************

********************************************************************************
*SET PATHS
********************************************************************************

cd "J:\EHR-Working\Krishnan\Kirsty\BC_respiratory_outcomes\paths"

do paths.do

********************************************************************************
*DATA MANAGEMENT 
********************************************************************************
//All events for all respiratory outcomes
cd $dofiles_respiratory
do cr_all_events_respiratory_all_codelists.do

//Combine neccesary files for creating exacerbation algorithms
cd $dofiles_respiratory
do cr_all_events_respiratory_outcomes_aurum.do

********************************************************************************
*DATA MANAGEMENT - LISTPAT FOR FIRST DX 
********************************************************************************
//Patient lists for first dx respiratory outcomes
cd $dofiles_respiratory
do cr_listpat_respiratory_firstdx_aurum_new.do

********************************************************************************
*DATA MANAGEMENT - LISTPAT FOR EXACERBATIONS 
********************************************************************************

//Patient lists chronic respiratory outcomes
cd $dofiles_respiratory

//Chronic diseases:
do cr_listpat_exacerbations.do // Main analaysis - 14 days washout period
do cr_listpat_exacerbations_30.do // For sensitivity analysis - 30 days washout period





