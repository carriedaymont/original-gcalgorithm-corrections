*Stata Code - Growth Data Cleaning 

/*
Code documentation is provided in Supplemental File #3. Numbered headings in this code correspond to 
numbered headings in that document. The following reference datasets needed to run the algorithm are available 
in Supplemental File #4: growthfile_cdc_ext, tanner_ht_vel_rev, who_ht_vel_3sd, who_ht_maxvel_3sd.

Please note that this Stata code currently runs only on Mac and does not work on Windows. 

*/

import growthdataforcleaning.csv, clear

************************************************************************************************
* #2: Data Setup
************************************************************************************************

*Drop duplicates for all variables
*817 dropped
duplicates drop

*Replace wt/ht as missing if ==0
replace wt=. if wt==0
replace ht=. if ht==0

*Add observation id variable for later matching
sort subjid agedays param measurement
gen obsid=string(_n, "%09.0f")
replace obsid = "I_" + obsid + "_I"
lab var obsid "Unique observation identifier"

save obsid, replace

*Create CDC z-scores and sd-scores

use obsid, clear
merge m:1 agedays sex using growthfile_cdc_ext
drop if _merge==2
drop _merge
drop cdc_hc* cdc_bmi_*

foreach p in wt ht{
	gen cdc`p'z=((`p'/cdc_`p'_m)^cdc_`p'_l-1)/(cdc_`p'_l*cdc_`p'_s)
	gen cdc`p'sd=(`p'-cdc_`p'_m)/cdc_`p'_csd_pos
	replace cdc`p'sd=(`p'-cdc_`p'_m)/cdc_`p'_csd_neg if `p'<cdc_`p'_m
}

lab var cdcwtz "CDC z-score for weight-for-age"
lab var cdchtz "CDC z-score for length/height-for-age"
lab var cdcwtsd "CDC SD score for weight-for-age"
lab var cdchtsd "CDC SD score for length/height-for-age"

save temp_sd, replace

************************************************************************************************
* #3: Re-centering
************************************************************************************************

*Use medians for years of age (combine years 19 and above)
*Combine sexes because very similar

use temp_sd, clear

gen year=int(agedays/365.25)
replace year=19 if year>=19
gen day_midyear=1 if agedays==int(year*365.25+365.25/2)
egen meanage=mean(agedays)

foreach p in wt ht {
	bysort param year: egen median`p'sd=median(cdc`p'sd)
	replace median`p'sd=. if day_midyear!=1
	ipolate median`p'sd agedays, gen(rc`p'sd)
	gen agedays_rc`p'=agedays if rc`p'sd!=.
	sort agedays_rc`p'
	gen rc`p'sd_minage=rc`p'sd[1]
	gsort -agedays_rc`p'
	gen rc`p'sd_maxage=rc`p'sd[1]
	replace rc`p'sd=rc`p'sd_minage if rc`p'sd==. & agedays<meanage
	replace rc`p'sd=rc`p'sd_maxage if rc`p'sd==. & agedays>meanage
	gen tbc`p'sd=cdc`p'sd-rc`p'sd
}



drop median*sd rc*sd_minage rc*sd_maxage meanage

lab var rcwtsd "Recentering value for weight SD"
lab var rchtsd "Recentering value for height SD"
lab var tbcwtsd "Weight SD score to be cleaned"
lab var tbchtsd "Height SD score to be cleaned"

*Set up exclusion variable and parameter-specific subject ID
foreach p in wt ht{
	gen exc_`p'=0
	replace exc_`p'=1 if `p'==.
}
lab var exc_wt "Exclusion marker for weight"
lab var exc_ht "Exclusion marker for height"

lab def exc 0 "Include" 1 "Missing"
lab val exc_* exc

foreach p in wt ht {
	gen str subjid_`p'=subjid
replace subjid_`p'="" if exc_`p'>0
}
lab var subjid_wt "Subject identifier for evaluable weights"
lab var subjid_ht "Subject identifier for evaluable heights"

order subjid

bysort subjid: gen tot=_N
save temp_rcsd_d, replace

************************************************************************************************
* #4: Dataset split
************************************************************************************************

use temp_rcsd_d, clear
keep if tot>30
drop tot
save temp_rcsd_d_long, replace
use temp_rcsd_d, clear
keep if tot>10 & tot<=30
drop tot
save temp_rcsd_d_medium, replace
use temp_rcsd_d, clear
keep if tot<=10
drop tot
save temp_rcsd_d_short, replace


************************************************************************************************
* #5: Temporary Duplicates
************************************************************************************************
*Temporarily select which value of duplicates to use based on median value for subject and parameter
*After removing larger errors will determine definitive value of duplicates to use

foreach length in short medium long {

use temp_rcsd_d_`length', clear

*Calculate medians for both parameters
foreach p in wt ht {
	duplicates tag subjid_`p' param agedays, gen(dup_`p')
	bysort subjid_`p': egen median_`p'sd_i=median(tbc`p'sd) if exc_`p'==0 & dup_`p'==0
	bysort subjid_`p': egen median_`p'sd=min(median_`p'sd_i)
}

foreach p in wt ht {
	if `p'==wt local o ht
	if `p'==ht local o wt
	gen absd_median_`p'sd=abs(tbc`p'sd-median_`p'sd)
	gen absd_median_`o'sd=abs(tbc`p'sd-median_`o'sd)
	sort subjid_`p' agedays absd_median_`p'sd absd_median_`o'sd
	by subjid_`p' agedays: gen d_sort_`p'_n=_n if exc_`p'==0
	replace exc_`p'=2 if exc_`p'==0 & d_sort_`p'_n>1
	replace subjid_`p'="" if exc_`p'>0
	drop absd_median_*  d_sort*
}

drop dup_* median*sd* 
lab def exc 2 "Duplicate exclusion (temporary)", add

save temp_dup_`length', replace
}

************************************************************************************************
* #7: Switches
************************************************************************************************


foreach length in short medium long {

use temp_dup_`length', clear

*Generate switched values for tbc*sd
foreach p in wt ht {
	gen `p'_orig=`p'
	if `p'==wt local o ht
	if `p'==ht local o wt
	gen `o'_i=`o' if exc_`o'==0
	bysort subjid agedays: egen `o'_b=min(`o'_i)
	gen tbc`p'sd_sw=(`o'_b-cdc_`p'_m)/cdc_`p'_csd_pos - rc`p'sd 
	replace tbc`p'sd_sw=((`o'_b-cdc_`p'_m)/cdc_`p'_csd_neg) - rc`p'sd if `o'_b<cdc_`p'_m 
}

*EWMA 
local e -1.5
local a 5
gen str3 ba=""
foreach p in wt ht {
	sort subjid_`p' agedays
	by subjid_`p': gen tot_`p'=_N
	by subjid_`p': gen sn_`p'=_n
	summ tot_`p' if exc_`p'==0
	local m=r(max)
	forvalues n=1/`m' {
		by subjid_`p': gen temp_tbc`p'sd_`n'=tbc`p'sd[`n'] 
		by subjid_`p': gen temp_age_`p'_`n'=agedays[`n'] 
		gen double temp_ewt_`p'_`n'=abs(`a'+agedays-temp_age_`p'_`n')^`e' if sn_`p'!=`n'
		gen double temp_ewtsd_`p'_`n'=temp_tbc`p'sd_`n'*temp_ewt_`p'_`n' 
		}
	egen double sumewt_`p'=rowtotal(temp_ewt_`p'_*) if exc_`p'==0
	egen double sumewtsd_`p'=rowtotal(temp_ewtsd_`p'_*) if exc_`p'==0
	gen ewma_`p'=sumewtsd_`p'/sumewt_`p' if subjid_`p'!=""
	gen dewma_`p'=tbc`p'sd-ewma_`p'
	foreach suf in sw {
		gen dewma_`p'_`suf'=tbc`p'sd_`suf'-ewma_`p'
	}
	*Before & After
	foreach t in bef aft {
		replace ba="`t'"
		gen double sumewt_`p'_`t'=sumewt_`p'
		gen double sumewtsd_`p'_`t'=sumewtsd_`p'
		gen ewma_`p'_`t'=ewma_`p'
		gen dewma_`p'_`t'=dewma_`p'
		foreach suf in sw {
		gen dewma_`p'_`suf'_`t'=dewma_`p'_`suf'
		}
		forvalues d=1/`m' {
			if ba=="bef" {
				local s=`d'+1
			}
			if ba=="aft" {
				local s=`d'-1
			}
			replace sumewt_`p'_`t'=sumewt_`p'-temp_ewt_`p'_`d' if temp_ewt_`p'_`d'!=. & sn_`p'==`s'
			replace sumewtsd_`p'_`t'=sumewtsd_`p'-temp_ewtsd_`p'_`d' if temp_ewtsd_`p'_`d'!=. & sn_`p'==`s'
			replace ewma_`p'_`t'=sumewtsd_`p'_`t'/sumewt_`p'_`t' if subjid_`p'!="" & sn_`p'==`s'
			replace dewma_`p'_`t'=tbc`p'sd-ewma_`p'_`t' if subjid_`p'!="" & sn_`p'==`s'
			foreach suf in sw {
				replace dewma_`p'_`suf'_`t'=tbc`p'sd_`suf'-ewma_`p'_`t' if subjid_`p'!="" & sn_`p'==`s'
			}
		}
	}	
	*Clean Up
	drop temp*
	sort subjid_`p' agedays
}

*Identify switches that meet criteria
gen switch_i=0
replace switch_i=1 if exc_wt==0 & tbcwtsd>4 & abs(tbcwtsd_sw)<3 & dewma_wt>3 & dewma_wt_bef>2 & dewma_wt_aft>2 & ///
	abs(dewma_wt_sw)<0.3 & abs(dewma_wt_sw_bef)<0.5 & abs(dewma_wt_sw_aft)<0.5
replace switch_i=1 if exc_ht==0 & tbchtsd<-7 & abs(tbchtsd_sw)<3 & dewma_ht<-6 & dewma_ht_bef<-5 & dewma_ht_aft<-5 & ///
	abs(dewma_ht_sw)<0.3 & abs(dewma_ht_sw_bef)<0.5 & abs(dewma_ht_sw_aft)<0.5
replace switch_i=0 if sn_wt==1|sn_wt==tot_wt
replace switch_i=0 if sn_ht==1|sn_ht==tot_ht

bysort subjid agedays: egen switch=total(switch_i) if exc_wt==0|exc_ht==0
recode switch 1=0 2=1
replace switch=0 if switch==.

*Rearrange switched values
foreach p in wt ht {
	if `p'==wt local o ht
	if `p'==ht local o wt
	replace tbc`p'sd=tbc`p'sd_sw if switch==1
	replace `p'=`o'_orig if switch==1
	replace tbc`p'sd=tbc`p'sd_sw if switch==1
}

replace exc_wt=1 if param=="WEIGHTKG" & switch==1
replace exc_ht=0 if param=="WEIGHTKG" & switch==1
replace exc_wt=0 if param=="HEIGHTCM" & switch==1
replace exc_ht=1 if param=="HEIGHTCM" & switch==1

foreach p in wt ht {
	replace subjid_`p'="" if exc_`p'==1
	replace subjid_`p'=subjid if exc_`p'==0
}

lab var wt_orig "Original weight value, no transformation"
lab var ht_orig "Original height value, no transformation"
lab var switch "Weight and height switched"
lab def switch 0 "No switch" 1 "Switch" 
lab val switch switch

drop ba-switch_i
drop *_i *_b *_sw

save temp_switch_`length', replace

}


************************************************************************************************
* #8: Unit errors
************************************************************************************************


foreach length in short medium long {

use temp_switch_`length', clear

*Create transformed values and tbc*sd scores
gen wt_d_2=wt/2.204622
gen wt_t_2=wt*2.204622
gen ht_d_2=ht/2.54
gen ht_t_2=ht*2.54

lab var wt_d_2 "Weight divided by 2.204622 to look for unit errors"
lab var wt_t_2 "Weight multiplied by 2.204622 to look for unit errors"
lab var ht_d_2 "Height divided by 2.54 to look for unit errors"
lab var ht_t_2 "Height multiplied by 2.54 to look for unit errors"


foreach p in wt ht {
	foreach dt in _d_2 _t_2 {
		gen tbc`p'sd`dt'=(`p'`dt'-cdc_`p'_m)/cdc_`p'_csd_pos - rc`p'sd
		replace tbc`p'sd`dt'=((`p'`dt'-cdc_`p'_m)/cdc_`p'_csd_neg)  - rc`p'sd if `p'`dt'<cdc_`p'_m
	}
}

*EWMA
local e -1.5
local a 5
gen str3 ba=""
foreach p in wt ht {
	sort subjid_`p' agedays
	by subjid_`p': gen tot_`p'=_N
	by subjid_`p': gen sn_`p'=_n
	summ tot_`p' if exc_`p'==0
	local m=r(max)
	forvalues n=1/`m' {
		by subjid_`p': gen temp_tbc`p'sd_`n'=tbc`p'sd[`n'] 
		by subjid_`p': gen temp_age_`p'_`n'=agedays[`n'] 
		gen double temp_ewt_`p'_`n'=abs(`a'+agedays-temp_age_`p'_`n')^`e' if sn_`p'!=`n'
		gen double temp_ewtsd_`p'_`n'=temp_tbc`p'sd_`n'*temp_ewt_`p'_`n' 
		}
	egen double sumewt_`p'=rowtotal(temp_ewt_`p'_*) if exc_`p'==0
	egen double sumewtsd_`p'=rowtotal(temp_ewtsd_`p'_*) if exc_`p'==0
	gen ewma_`p'=sumewtsd_`p'/sumewt_`p' if subjid_`p'!=""
	gen dewma_`p'=tbc`p'sd-ewma_`p'
	foreach dt in _d_2 _t_2 {
		gen dewma_`p'`dt'=tbc`p'sd`dt'-ewma_`p'
	}
	*Before & After
	foreach t in bef aft {
		replace ba="`t'"
		gen double sumewt_`p'_`t'=sumewt_`p'
		gen double sumewtsd_`p'_`t'=sumewtsd_`p'
		gen ewma_`p'_`t'=ewma_`p'
		gen dewma_`p'_`t'=dewma_`p'
		foreach dt in _d_2 _t_2 {
			gen dewma_`p'`dt'_`t'=dewma_`p'`dt'
		}
		forvalues d=1/`m' {
			if ba=="bef" {
				local s=`d'+1
			}
			if ba=="aft" {
				local s=`d'-1
			}
			replace sumewt_`p'_`t'=sumewt_`p'-temp_ewt_`p'_`d' if temp_ewt_`p'_`d'!=. & sn_`p'==`s'
			replace sumewtsd_`p'_`t'=sumewtsd_`p'-temp_ewtsd_`p'_`d' if temp_ewtsd_`p'_`d'!=. & sn_`p'==`s'
			replace ewma_`p'_`t'=sumewtsd_`p'_`t'/sumewt_`p'_`t' if subjid_`p'!="" & sn_`p'==`s'
			replace dewma_`p'_`t'=tbc`p'sd-ewma_`p'_`t' if subjid_`p'!="" & sn_`p'==`s'
			foreach dt in _d_2 _t_2 {
				replace dewma_`p'`dt'_`t'=tbc`p'sd`dt'-ewma_`p'_`t' if subjid_`p'!="" & sn_`p'==`s'
			}
		}
	}	
	*Clean Up
	drop temp*
	sort subjid_`p' agedays
	by subjid_`p': gen d_nextsd_`p'=tbc`p'sd-tbc`p'sd[_n+1]
	by subjid_`p': gen d_prevsd_`p'=tbc`p'sd-tbc`p'sd[_n-1]
}

*Identify unit errors that meet criteria
gen replaced_wt_d_2=1 if dewma_wt>3 & dewma_wt_bef>2 & dewma_wt_aft>2 & tbcwtsd>3 & ///
	d_nextsd_wt>2 & d_nextsd_wt!=. & d_prevsd_wt>2 & d_prevsd_wt!=. &   ///
	abs(dewma_wt_d_2)<0.3 & abs(dewma_wt_d_2_bef)<0.5 & abs(dewma_wt_d_2_aft)<0.5 & abs(tbcwtsd_d_2)<3 & exc_wt==0
gen replaced_wt_t_2=1 if dewma_wt<-3 & dewma_wt_bef<-2 & dewma_wt_aft<-2 & tbcwtsd<-3 & ///
	d_nextsd_wt<-2 & d_nextsd_wt!=. & d_prevsd_wt<-2 & d_prevsd_wt!=. &   ///
	abs(dewma_wt_t_2)<0.3 & abs(dewma_wt_t_2_bef)<0.5 & abs(dewma_wt_t_2_aft)<0.5 & abs(tbcwtsd_t_2)<3 & exc_wt==0	
gen replaced_ht_d_2=1 if dewma_ht>5 & dewma_ht_bef>4 & dewma_ht_aft>4 & tbchtsd>7 & ///
	d_nextsd_ht>4 & d_nextsd_ht!=. & d_prevsd_ht>4 & d_prevsd_ht!=. &   ///
	abs(dewma_ht_d_2)<0.3 & abs(dewma_ht_d_2_bef)<0.5 & abs(dewma_ht_d_2_aft)<0.5 & abs(tbchtsd_d_2)<3 & exc_ht==0
gen replaced_ht_t_2=1 if dewma_ht<-5 & dewma_ht_bef<-4 & dewma_ht_aft<-4 & tbchtsd<-7 & ///
	d_nextsd_ht<-4 & d_nextsd_ht!=. & d_prevsd_ht<-4 & d_prevsd_ht!=. &   ///
	abs(dewma_ht_t_2)<0.3 & abs(dewma_ht_t_2_bef)<0.5 & abs(dewma_ht_t_2_aft)<0.5 & abs(tbchtsd_t_2)<3 & exc_ht==0	

*Adjust wt/ht and tbc*sd
foreach p in wt ht {
	foreach dt in _d_2 _t_2 {
		replace tbc`p'sd=tbc`p'sd`dt' if replaced_`p'`dt'==1
		replace `p'=`p'`dt' if replaced_`p'`dt'==1
	}
}
gen unit_transformed=0
replace unit_transformed=1 if replaced_wt_d_2==1|replaced_wt_t_2==1|replaced_ht_d_2==1|replaced_ht_t_2==1

lab var tbcwtsd_d_2 "SD score for weight transformed by dividng by 2.204622"
lab var tbcwtsd_t_2 "SD score for weight transformed by multiplying by 2.204622"
lab var tbchtsd_d_2 "SD score for height transformed by dividng by 2.54"
lab var tbchtsd_t_2 "SD score for height transformed by multiplying by 2.54"
drop ba-d_prevsd_ht
lab var replaced_wt_d_2 "If 1, wt & SD weight score replaced with transformed value (wt_orig/2.204622)"
lab var replaced_wt_t_2 "If 1, wt & SD weight score replaced with transformed value (wt_orig*2.204622)"
lab var replaced_ht_d_2 "If 1, ht &  SD height score replaced with transformed value (ht_orig/2.54)"
lab var replaced_ht_t_2 "If 1, ht &  SD height score replaced with transformed value (ht_orig*2.54)"
lab var unit_transformed "If 1, likely unit error, value and SD score replaced with transformed value"
lab def replaced 0 "Original value kept" 1 "Transformed value used"
lab val replaced_wt_d_2 replaced
lab val replaced_wt_t_2 replaced
lab val replaced_ht_d_2 replaced
lab val replaced_ht_t_2 replaced
lab val unit_transformed replaced
lab var wt_orig "Original weight value, no transformation"
lab var ht_orig "Original height value, no transformation"




save temp_unit_`length', replace

}


************************************************************************************************
* #9: Carried forward
************************************************************************************************


foreach length in  short medium long  {

use temp_unit_`length', clear

foreach p in wt ht {
	*Identify carried forward if no duplicates
	sort subjid_`p' agedays
	by subjid_`p': gen d_prev_`p'=measurement-measurement[_n-1]
	replace exc_`p'=3 if d_prev_`p'==0 & exc_`p'!=1 & unit_transformed==0 & switch==0
	replace subjid_`p'="" if exc_`p'>0
	drop d_prev_`p'
	*Identify carried forward if duplicates
	gen subjid_`p'2=subjid if exc_`p'==0|exc_`p'==2|exc_`p'==3
	sort subjid_`p'2 agedays exc_`p'
	by subjid_`p'2 agedays: gen tempn=_n if subjid_`p'2!=""
	by subjid_`p'2 agedays: gen temptot=_N if subjid_`p'2!=""
	sort subjid_`p'2 tempn agedays
	by subjid_`p'2 tempn: gen agen_i=_n if tempn==1 & subjid_`p'2!=""
	bysort subjid_`p'2 agedays: egen agen=min(agen_i)
	drop agen_i
	gen agen_prev=agen-1 if agen>1
	sort subjid_`p'2 tempn agen
	by subjid_`p'2 tempn: gen tempn_prev_i=tempn[_n-1]
	bysort subjid_`p'2 agen: egen tempn_prev=min(tempn_prev_i)
	drop tempn_prev_i
	egen maxtempn=max(tempn)
	egen maxagen=max(agen)
	local mtempn=maxtempn[1]
	local magen=maxagen[1]
	forvalues tn=1/`mtempn' {
		gen d_prevn_`tn'=.
	}
	forvalues an=2/`magen' {
		forvalues tn=1/`mtempn' {
			gen tempcomp_i=measurement if agen==(`an'-1) & tempn==`tn' & subjid_`p'2!=""
			bysort subjid_`p'2: egen tempcomp=min(tempcomp_i)
			replace d_prevn_`tn'=abs(measurement-tempcomp) if agen==`an' & subjid_`p'2!=""
			drop tempcomp tempcomp_i
		}
	}
	egen mincomp=rowmin(d_prevn_1-d_prevn_`mtempn') if subjid_`p'2!="" & agen>=2
	replace exc_`p'=3 if mincomp==0 & exc_`p'!=1 & unit_transformed==0 & switch==0
	replace subjid_`p'="" if exc_`p'>0
	drop subjid_`p'2-mincomp
}

*Redo temporary duplicates
foreach p in wt ht {
	duplicates tag subjid_`p' param agedays, gen(dup_`p')
	bysort subjid_`p': egen median_`p'sd_i=median(tbc`p'sd) if exc_`p'==0 & dup_`p'==0
	bysort subjid_`p': egen median_`p'sd=min(median_`p'sd_i)
}

foreach p in wt ht {
	if `p'==wt local o ht
	if `p'==ht local o wt
	gen absd_median_`p'sd=abs(tbc`p'sd-median_`p'sd)
	gen absd_median_`o'sd=abs(tbc`p'sd-median_`o'sd)
	sort subjid_`p' agedays absd_median_`p'sd absd_median_`o'sd
	by subjid_`p' agedays: gen d_sort_`p'_n=_n if exc_`p'==0
	replace exc_`p'=2 if exc_`p'==0 & d_sort_`p'_n>1
	replace subjid_`p'="" if exc_`p'>0
	drop absd_median_*  d_sort*
}

drop dup_* median_*
lab def exc 3 "Carried forward", add 

save temp_forward_`length', replace

}


************************************************************************************************
* #10: Extreme errors with SD cutoffs
************************************************************************************************


foreach length in short medium long  {

use temp_forward_`length', clear

*Identify extreme errors with SD and Z-score cutoffs
foreach p in wt ht {
	replace exc_`p'=4 if abs(tbc`p'sd)>25 & (exc_`p'==0|exc_`p'==2)
	replace exc_`p'=4 if abs(cdc`p'z)>25 & switch==0 & unit_transformed==0 & (exc_`p'==0|exc_`p'==2)
	replace subjid_`p'="" if exc_`p'>0
}

*Redo temporary duplicates
foreach p in wt ht {
	duplicates tag subjid_`p' param agedays, gen(dup_`p')
	bysort subjid_`p': egen median_`p'sd_i=median(tbc`p'sd) if exc_`p'==0 & dup_`p'==0
	bysort subjid_`p': egen median_`p'sd=min(median_`p'sd_i)
}


foreach p in wt ht {
	if `p'==wt local o ht
	if `p'==ht local o wt
	gen absd_median_`p'sd=abs(tbc`p'sd-median_`p'sd)
	gen absd_median_`o'sd=abs(tbc`p'sd-median_`o'sd)
	sort subjid_`p' agedays absd_median_`p'sd absd_median_`o'sd
	by subjid_`p' agedays: gen d_sort_`p'_n=_n if exc_`p'==0
	replace exc_`p'=2 if exc_`p'==0 & d_sort_`p'_n>1
	replace subjid_`p'="" if exc_`p'>0
	drop absd_median_*  d_sort*
}

drop dup* median_* 


************************************************************************************************
* #11: Extreme errors with EWMA
************************************************************************************************

*EWMA
foreach p in wt ht {	
	local e -1.5
	local a 5
	local i 1 
	while `i'>0 {
		gen str3 ba=""
		gen seldup_`p'=1 if exc_`p'==0
		replace exc_`p'=0 if exc_`p'==2
		replace subjid_`p'=subjid if exc_`p'==0
		sort subjid_`p' agedays
		by subjid_`p': gen tot_`p'=_N
		by subjid_`p': gen sn_`p'=_n
		summ tot_`p' if exc_`p'==0
		local m=r(max)
		forvalues n=1/`m' {
			by subjid_`p': gen temp_tbc`p'sd_`n'=tbc`p'sd[`n'] if seldup_`p'[`n']==1
			by subjid_`p': gen temp_age_`p'_`n'=agedays[`n'] if seldup_`p'[`n']==1
			gen double temp_ewt_`p'_`n'=abs(`a'+agedays-temp_age_`p'_`n')^`e' if (agedays-temp_age_`p'_`n')!=0 
			gen double temp_ewtsd_`p'_`n'=temp_tbc`p'sd_`n'*temp_ewt_`p'_`n' 
			}
		egen double sumewt_`p'=rowtotal(temp_ewt_`p'_*) if exc_`p'==0
		egen double sumewtsd_`p'=rowtotal(temp_ewtsd_`p'_*) if exc_`p'==0
		gen ewma_`p'=sumewtsd_`p'/sumewt_`p' if subjid_`p'!=""
		gen dewma_`p'=tbc`p'sd-ewma_`p'
		*Before & After
		foreach t in bef aft {
			replace ba="`t'"
			gen double sumewt_`p'_`t'=sumewt_`p'
			gen double sumewtsd_`p'_`t'=sumewtsd_`p'
			gen ewma_`p'_`t'=ewma_`p'
			gen dewma_`p'_`t'=dewma_`p'
			forvalues d=1/`m' {
				if ba=="bef" {
					local s=`d'+1
				}
				if ba=="aft" {
					local s=`d'-1
				}
				replace sumewt_`p'_`t'=sumewt_`p'-temp_ewt_`p'_`d' if temp_ewt_`p'_`d'!=. & sn_`p'==`s'
				replace sumewtsd_`p'_`t'=sumewtsd_`p'-temp_ewtsd_`p'_`d' if temp_ewtsd_`p'_`d'!=. & sn_`p'==`s'
				replace ewma_`p'_`t'=sumewtsd_`p'_`t'/sumewt_`p'_`t' if subjid_`p'!="" & sn_`p'==`s'
				replace dewma_`p'_`t'=tbc`p'sd-ewma_`p'_`t' if subjid_`p'!="" & sn_`p'==`s'
			}
		}	
		*Clean Up
		drop temp*
		sort subjid_`p' agedays
		by subjid_`p': gen d_nextsd_`p'=tbc`p'sd-tbc`p'sd[_n+1]
		by subjid_`p': gen d_prevsd_`p'=tbc`p'sd-tbc`p'sd[_n-1]
		gen abs_tbc`p'sd=abs(tbc`p'sd)
		bysort subjid_`p': egen max_abs_tbc`p'sd=max(abs_tbc`p'sd)
		gen abssum_sd_dewma_`p'=abs(tbc`p'sd+dewma_`p')
		gen extreme_exc_`p'=1 if dewma_`p'>3.5 & (dewma_`p'_bef>3|dewma_`p'==.) & (dewma_`p'_aft>3|dewma_`p'_aft==.) & tbc`p'sd>3.5 & exc_`p'==0 & tot_`p'>2
		replace extreme_exc_`p'=1 if dewma_`p'<-3.5 & (dewma_`p'_bef<-3|dewma_`p'==.) & (dewma_`p'_aft<-3|dewma_`p'_aft==.) & tbc`p'sd<-3.5 & exc_`p'==0 & tot_`p'>2
		gsort subjid_`p' extreme_exc_`p' -abssum_sd_dewma_`p'
		by subjid_`p' extreme_exc_`p': gen ext_exc_num_`p'=_n if extreme_exc_`p'==1
		replace ext_exc_num_`p'=0 if ext_exc_num_`p'==.
		egen max_ext_exc_num_`p'=max(ext_exc_num_`p')
		replace exc_`p'=5 if ext_exc_num_`p'==1
		replace exc_`p'=6 if dewma_`p'>3.5 & tbc`p'sd>3.5 & max_abs_tbc`p'sd==abs(tbc`p'sd) & exc_`p'==0 & tot_`p'==2 
		replace exc_`p'=6 if dewma_`p'<-3.5 & tbc`p'sd<-3.5 & max_abs_tbc`p'sd==abs(tbc`p'sd) & exc_`p'==0 & tot_`p'==2 
		replace subjid_`p'="" if exc_`p'>0
		local i=max_ext_exc_num_`p'[1]
		drop ba-max_ext_exc_num_`p'
		*Redo duplicates
		foreach psub in wt ht {
			duplicates tag subjid_`psub' param agedays, gen(dup_`psub')
			bysort subjid_`psub': egen median_`psub'sd_i=median(tbc`psub'sd) if exc_`psub'==0 & dup_`psub'==0
			bysort subjid_`psub': egen median_`psub'sd=min(median_`psub'sd_i)
		}
		if `p'==wt local o ht
		if `p'==ht local o wt
		gen absd_median_`p'sd=abs(tbc`p'sd-median_`p'sd)
		gen absd_median_`o'sd=abs(tbc`p'sd-median_`o'sd)
		sort subjid_`p' agedays absd_median_`p'sd absd_median_`o'sd
		by subjid_`p' agedays: gen d_sort_`p'_n=_n if exc_`p'==0
		replace exc_`p'=2 if exc_`p'==0 & d_sort_`p'_n>1
		replace subjid_`p'="" if exc_`p'>0
		drop dup* median_* absd_* d_sort*
	}
}

lab def exc 4 "|SD|>25" 5 "EWMA1 >=3" 6 "EWMA1 2", add

save temp_ewma1_`length', replace
}

************************************************************************************************
* #12: Duplicates with EWMA
************************************************************************************************


foreach length in short medium long  {

use temp_ewma1_`length', clear

foreach p in wt ht {
	replace exc_`p'=0 if exc_`p'==2
	replace subjid_`p'=subjid if exc_`p'==0
	duplicates tag subjid_`p' agedays, gen(dup_`p')
	bysort subjid_`p': egen median_`p'sd_i=median(tbc`p'sd) if exc_`p'==0 & dup_`p'==0
	bysort subjid_`p': egen median_`p'sd=min(median_`p'sd)	
}

recode wt (0/9.9999=0.25) (10/29.99999=0.5) (30/1000=1), gen(wtmult_i)
bysort subjid param agedays: egen wtmult=min(wtmult_i)
drop wtmult_i


*EWMA
local e -1.5
local a 5
gen str3 ba=""
foreach p in wt ht {
	if `p'==wt local o ht
	if `p'==ht local o wt
	replace median_`p'sd=0 if median_`o'sd==.
	gen absd_median_`p'sd=abs(tbc`p'sd-median_`p'sd)
	gen absd_median_`o'sd=abs(tbc`p'sd-median_`o'sd)
	sort subjid_`p' agedays absd_median_`p'sd absd_median_`o'sd
	by subjid_`p' agedays: gen `p'_sort_n=_n if exc_`p'==0
	by subjid_`p': gen tot_`p'=_N
	by subjid_`p': gen sn_`p'=_n
	summ tot_`p' if exc_`p'==0
	local m=r(max)
	forvalues n=1/`m' {
		by subjid_`p': gen temp_tbc`p'sd_`n'=tbc`p'sd[`n'] if `p'_sort_n[`n']==1
		by subjid_`p': gen temp_age_`p'_`n'=agedays[`n'] if `p'_sort_n[`n']==1
		gen double temp_ewt_`p'_`n'=abs(`a'+agedays-temp_age_`p'_`n')^`e' if (agedays-temp_age_`p'_`n')!=0 
		gen double temp_ewtsd_`p'_`n'=temp_tbc`p'sd_`n'*temp_ewt_`p'_`n' 
		}
	egen double sumewt_`p'=rowtotal(temp_ewt_`p'_*) 
	egen double sumewtsd_`p'=rowtotal(temp_ewtsd_`p'_*)
	gen ewma_`p'=sumewtsd_`p'/sumewt_`p' if subjid_`p'!=""
	gen dewma_`p'=tbc`p'sd-ewma_`p'	
	*Clean Up
	gen abssum_`p'=abs(2*dewma_`p')+abs(tbc`p'sd)
	replace abssum_`p'=abs(tbc`p'sd) if dewma_`p'==.
	sort subjid_`p' agedays abssum_`p'
	by subjid_`p' agedays: gen `p'_ewma_n=_n if exc_`p'==0 
	replace exc_`p'=7 if `p'_ewma_n>1 & exc_`p'==0 
	*Exclude all duplicates if too difficult to tell which is likely representative
	sort subjid_`p' agedays `p'
	by subjid_`p' agedays: gen d`p'=`p'[dup_`p'+1]-`p'[1] if dup_`p'>0
	by subjid_`p' agedays: gen dup_`p'_1=1 if _n==2
	by subjid_`p': egen totdup_`p'=total(dup_`p'_1)
	gen nodup_`p'=dup_`p'==0 if exc_`p'==0
	bysort subjid_`p': egen totnodup_`p'=total(nodup_`p')
	if `p'==ht local d 3
	if `p' ==wt local d wtmult
	gen propdup_`p'=(totdup_`p'/(totdup_`p'+nodup_`p'))
	replace exc_`p'=7 if exc_`p'==0 & dup_`p'>0 & (totdup_`p'/(totdup_`p'+totnodup_`p'))>0.5 & totnodup_`p'>0 & d`p'>`d'
	replace exc_`p'=7 if exc_`p'==0 & dup_`p'>0 & totnodup_`p'==0 & d`p'>`d'
	replace subjid_`p'="" if exc_`p'>0
	drop absd_median* 
}

*Replace exc_*=7 if previously excluded as extremes
foreach p in wt ht {
	gen subjid_`p'2=subjid_`p'
	replace subjid_`p'2=subjid if exc_`p'>=4 & exc_`p'<=6
	duplicates tag subjid_`p'2 agedays, gen(dup_`p'2)
	replace exc_`p'=7 if exc_`p'>=4 & exc_`p'<=6 & dup_`p'2>0
	replace subjid_`p'="" if exc_`p'>0
	drop subjid_`p'2 dup_`p'2
}


lab def exc 7 "Duplicate", modify
drop ba-totdup_ht
drop median* dup_wt dup_ht
drop wtmult



************************************************************************************************
* #13: Plus/minus measurements
************************************************************************************************


gen wt_plus=wt+0.05*wt
gen wt_minus=wt-0.05*wt
gen ht_plus=ht+1
gen ht_minus=ht-1

lab var wt_plus "Weight  +  5%"
lab var wt_minus "Weight - 5%"
lab var ht_plus "Height + 1 cm"
lab var ht_minus "Height - 1 cm"

foreach p in wt ht {
	foreach d in plus minus {
		gen tbc`p'sd_`d'=(`p'_`d'-cdc_`p'_m)/cdc_`p'_csd_pos - rc`p'sd
		replace tbc`p'sd_`d'=(`p'_`d'-cdc_`p'_m)/cdc_`p'_csd_neg if `p'_`d'<cdc_`p'_m - rc`p'sd
	}
}

lab var tbcwtsd_plus "SD score for weight + 5%"
lab var tbcwtsd_minus "SD score for weight + 5%"
lab var tbchtsd_plus "SD score for height + 1 cm"
lab var tbchtsd_minus "SD score for height - 1 cm"


save temp_dup2_`length', replace

}


************************************************************************************************
* #14: Moderate errors with EWMA
************************************************************************************************


foreach length in short medium long  {

use temp_dup2_`length', clear

*EWMA
local e -1.5
local a 5
foreach p in wt ht {
	local i 1
	while `i'>0 {
		gen str3 ba=""
		sort subjid_`p' agedays
		by subjid_`p': gen tot_`p'=_N
		by subjid_`p': gen sn_`p'=_n
		summ tot_`p' if exc_`p'==0
		local m=r(max)
		forvalues n=1/`m' {
			by subjid_`p': gen temp_tbc`p'sd_`n'=tbc`p'sd[`n'] 
			by subjid_`p': gen temp_age_`p'_`n'=agedays[`n'] 
			gen double temp_ewt_`p'_`n'=abs(`a'+agedays-temp_age_`p'_`n')^`e' if sn_`p'!=`n'
			gen double temp_ewtsd_`p'_`n'=temp_tbc`p'sd_`n'*temp_ewt_`p'_`n' 
			}
		egen double sumewt_`p'=rowtotal(temp_ewt_`p'_*) if exc_`p'==0
		egen double sumewtsd_`p'=rowtotal(temp_ewtsd_`p'_*) if exc_`p'==0
		gen ewma_`p'=sumewtsd_`p'/sumewt_`p' if subjid_`p'!=""
		gen dewma_`p'=tbc`p'sd-ewma_`p'
		*Before & After
		foreach t in bef aft {
			replace ba="`t'"
			gen double sumewt_`p'_`t'=sumewt_`p'
			gen double sumewtsd_`p'_`t'=sumewtsd_`p'
			gen ewma_`p'_`t'=ewma_`p'
			gen dewma_`p'_`t'=dewma_`p'
			forvalues d=1/`m' {
				if ba=="bef" {
					local s=`d'+1
				}
				if ba=="aft" {
					local s=`d'-1
				}
				replace sumewt_`p'_`t'=sumewt_`p'-temp_ewt_`p'_`d' if temp_ewt_`p'_`d'!=. & sn_`p'==`s'
				replace sumewtsd_`p'_`t'=sumewtsd_`p'-temp_ewtsd_`p'_`d' if temp_ewtsd_`p'_`d'!=. & sn_`p'==`s'
				replace ewma_`p'_`t'=sumewtsd_`p'_`t'/sumewt_`p'_`t' if subjid_`p'!="" & sn_`p'==`s'
				replace dewma_`p'_`t'=tbc`p'sd-ewma_`p'_`t' if subjid_`p'!="" & sn_`p'==`s'
			}
		}	
		drop temp*
		drop sn_`p'
		sort subjid_`p' agedays
		by subjid_`p': gen sn_`p'=_n
		by subjid_`p': gen vistot_`p'=_N
		by subjid_`p': gen d_agedays_prev_`p'=agedays-agedays[_n-1]
		by subjid_`p': gen d_agedays_next_`p'=agedays[_n+1]-agedays
		gen fl_`p'=0
		replace fl_`p'=1 if sn_`p'==1
		replace fl_`p'=2 if sn_`p'==vistot_`p' * (d_agedays_prev_`p'<(365.25*2))
		replace fl_`p'=3 if sn_`p'==vistot_`p' * (d_agedays_prev_`p'>=(365.25*2))
		replace fl_`p'=. if exc_`p'>0
		replace fl_`p'=. if vistot_`p'<3
		sort subjid_`p' agedays
		by subjid_`p': gen tbc`p'sd_prev=tbc`p'sd[_n-1]
		by subjid_`p': gen tbc`p'sd_next=tbc`p'sd[_n+1]
		by subjid_`p': gen d_prevsd_`p'=tbc`p'sd-tbc`p'sd[_n-1]
		by subjid_`p': gen d_nextsd_`p'=tbc`p'sd-tbc`p'sd[_n+1]
		foreach d in plus minus {
			by subjid_`p': gen d_prevsd_`d'_`p'=tbc`p'sd_`d'-tbc`p'sd[_n-1]
			by subjid_`p': gen d_nextsd_`d'_`p'=tbc`p'sd_`d'-tbc`p'sd[_n+1]
		}
		gen exc_temp_`p'=.
		*Identify possible exclusions
		replace exc_temp_`p'=8 if fl_`p'==0 & dewma_`p'>1 & dewma_`p'_bef>1 & dewma_`p'_aft>1 & d_nextsd_`p'>1 & d_prevsd_`p'>1 & ///
			d_prevsd_plus_`p'>1 & d_prevsd_minus_`p'>1 & d_nextsd_plus_`p'>1 & d_nextsd_minus_`p'>1
		replace exc_temp_`p'=8 if fl_`p'==0 & dewma_`p'<-1 & dewma_`p'_bef<-1 & dewma_`p'_aft<-1 & d_nextsd_`p'<-1 & d_prevsd_`p'<-1 & ///
			d_prevsd_plus_`p'<-1 & d_prevsd_minus_`p'<-1 & d_nextsd_plus_`p'<-1 & d_nextsd_minus_`p'<-1
		replace exc_temp_`p'=9 if fl_`p'==1 & d_agedays_next_`p'<365.25 & dewma_`p'>2 & dewma_`p'_aft>1 & d_nextsd_`p'>1 & ///
			d_nextsd_plus_`p'>1 & d_nextsd_minus_`p'>1
		replace exc_temp_`p'=9 if fl_`p'==1 & d_agedays_next_`p'<365.25 & dewma_`p'<-2 & dewma_`p'_aft<-1 & d_nextsd_`p'<-1 &  ///
			d_nextsd_plus_`p'<-1 & d_nextsd_minus_`p'<-1	
		replace exc_temp_`p'=10 if fl_`p'==1 & d_agedays_next_`p'>365.25 & dewma_`p'>3 & dewma_`p'_aft>1 & d_nextsd_`p'>1 & ///
			d_nextsd_plus_`p'>1 & d_nextsd_minus_`p'>1
		replace exc_temp_`p'=10 if fl_`p'==1 & d_agedays_next_`p'>365.25 & dewma_`p'<-3 & dewma_`p'_aft<-1 & d_nextsd_`p'<-1 & ///
			d_nextsd_plus_`p'<-1 & d_nextsd_minus_`p'<-1	
		replace exc_temp_`p'=11 if fl_`p'==2 & abs(tbc`p'sd_prev)<2 & dewma_`p'>2 & dewma_`p'_bef>1 & d_prevsd_`p'>1 & ///
			d_prevsd_plus_`p'>1 & d_prevsd_minus_`p'>1 
		replace exc_temp_`p'=11 if fl_`p'==2 & abs(tbc`p'sd_prev)<2 & dewma_`p'<-2 & dewma_`p'_bef<-1 & d_prevsd_`p'<-1 & ///
			d_prevsd_plus_`p'<-1 & d_prevsd_minus_`p'<-1 
		replace exc_temp_`p'=12 if fl_`p'==2 & abs(tbc`p'sd_prev)>=2 & dewma_`p'>abs(tbc`p'sd_prev) & dewma_`p'_bef>1 & d_prevsd_`p'>1 & ///
			d_prevsd_plus_`p'>1 & d_prevsd_minus_`p'>1 
		replace exc_temp_`p'=12 if fl_`p'==2 & abs(tbc`p'sd_prev)>=2 & dewma_`p'<(-1 *abs(tbc`p'sd_prev)) & dewma_`p'_bef<-1 & d_prevsd_`p'<-1 & ///
			d_prevsd_plus_`p'<-1 & d_prevsd_minus_`p'<-1  
		if `p'==wt local o ht
		if `p'==ht local o wt
		bysort subjid_`o': egen median_tbc`o'sd_i=median(tbc`o'sd) if exc_`o'==0
		bysort subjid: egen median_tbc`o'sd=min(median_tbc`o'sd_i)	
		gen comp_tbc`o'sd_i=tbc`o'sd if exc_`o'==0
		bysort subjid agedays: egen comp_tbc`o'sd=min(comp_tbc`o'sd_i)
		replace exc_temp_`p'=13 if fl_`p'==3 & abs(tbc`p'sd_prev)<2 & dewma_`p'>3 & dewma_`p'_bef>1 & d_prevsd_`p'>1 & ///
			d_prevsd_plus_`p'>1 & d_prevsd_minus_`p'>1 & ///
			(((tbc`p'sd-tbc`o'sd)>4 & tbc`o'sd!=.)|((tbc`p'sd-median_tbc`o'sd)>4 & median_tbc`o'sd!=. & tbc`o'sd==.)|median_tbc`o'sd==.)
		replace exc_temp_`p'=13 if fl_`p'==3 & abs(tbc`p'sd_prev)<2 & dewma_`p'<-3 & dewma_`p'_bef<-1 & d_prevsd_`p'<-1 & ///
			d_prevsd_plus_`p'<-1 & d_prevsd_minus_`p'<-1 & ///
			(((tbc`p'sd-tbc`o'sd)<-4 & tbc`o'sd!=.)|((tbc`p'sd-median_tbc`o'sd)<-4 & median_tbc`o'sd!=. & tbc`o'sd==.)|median_tbc`o'sd==.)
		replace exc_temp_`p'=14 if fl_`p'==3 & abs(tbc`p'sd_prev)>=2 & dewma_`p'>(1+abs(tbc`p'sd_prev)) & dewma_`p'_bef>1 & d_prevsd_`p'>1 & ///
			d_prevsd_plus_`p'>1 & d_prevsd_minus_`p'>1 & ///
			(((tbc`p'sd-tbc`o'sd)>4 & tbc`o'sd!=.)|((tbc`p'sd-median_tbc`o'sd)>4 & median_tbc`o'sd!=. & tbc`o'sd==.)|median_tbc`o'sd==.)
		replace exc_temp_`p'=14 if fl_`p'==3 & abs(tbc`p'sd_prev)>=2 & dewma_`p'<(-1-abs(tbc`p'sd_prev)) & dewma_`p'_bef<-1 & d_prevsd_`p'<-1 & ///
			d_prevsd_plus_`p'<-1 & d_prevsd_minus_`p'<-1 & ///
			(((tbc`p'sd-tbc`o'sd)<-4 & tbc`o'sd!=.)|((tbc`p'sd-median_tbc`o'sd)<-4 & median_tbc`o'sd!=. & tbc`o'sd==.)|median_tbc`o'sd==.)
		gen abssum_sd_dewma_`p'=abs(tbc`p'sd+dewma_`p')
		gen exc_ind_`p'=exc_temp_`p'!=.
		gsort subjid_`p' agedays exc_ind_`p' fl_`p' -abssum_sd_dewma_`p'
		by subjid_`p' agedays exc_ind_`p': gen exc_num_`p'=_n if exc_ind_`p'==1
		replace exc_num_`p'=0 if exc_num_`p'==.
		egen max_exc_num_`p'=max(exc_num_`p')
		*Identify exclusions for this round
		replace exc_`p'=exc_temp_`p' if exc_num_`p'==1
		replace subjid_`p'="" if exc_`p'>0
		local i=max_exc_num_`p'[1]
		drop ba-max_exc_num_`p'
	}
}
lab def exc 8 "EWMA2 middle" 9 "EWMA2 1st" 10 "EWMA2 1st ext" 11 "EWMA2 last" 12 "EWMA last ext" ///
	13 "EWMA2 last long" 14 "EWMA2 last long ext", add

save temp_ewma2_`length', replace
}

************************************************************************************************
* #15: Absolute differences
************************************************************************************************


foreach length in short medium long  {

use temp_ewma2_`length', clear

local e -1.5
local a 5
local p ht
local i 1 
while `i'>0 {
	gen str3 ba=""
	*Identify midpoints and intervals
	sort subjid_`p' agedays
	by subjid_`p': gen d_agedays_`p'=agedays[_n+1]-agedays if subjid_`p'!=""
	by subjid_`p': gen mid_agedays_`p'=0.5*(agedays+agedays[_n+1]) if subjid_`p'!=""
	gen tanner_months=6+12*(round(mid_agedays_ht/365.25))
	*Tanner values
	merge m:1 sex tanner_months using tanner_`p'_vel_rev.dta, nogen
	drop if subjid==""
	*Scale min and max for interval and add extra allowance
	gen mindiff_`p'=0.5*min_`p'_vel*(d_agedays/365.25)^2-3 if d_agedays<365.25
	replace mindiff_`p'=0.5*min_`p'_vel-3 if d_agedays>365.25
	gen maxdiff_`p'=2*max_`p'_vel*(d_agedays/365.25)^1.5+5.5 if d_agedays>365.25
	replace maxdiff_`p'=2*max_`p'_vel*(d_agedays/365.25)^0.33+5.5 if d_agedays<365.25
	*WHO values
	gen whoagegrp_`p'=round(agedays/30.4375, 1) if agedays<=(2*365.25) 
	gen whoinc_age_`p'=.
	replace whoinc_age_`p'=1 if d_agedays_`p'>=20 & d_agedays_`p'<46
	replace whoinc_age_`p'=2 if d_agedays_`p'>=46 & d_agedays_`p'<76
	replace whoinc_age_`p'=3 if d_agedays_`p'>=76 & d_agedays_`p'<107
	replace whoinc_age_`p'=4 if d_agedays_`p'>=107 & d_agedays_`p'<153
	replace whoinc_age_`p'=6 if d_agedays_`p'>=153 & d_agedays_`p'<199	
	merge m:1 sex whoagegrp_`p' using who_`p'_vel_3sd.dta, nogen
	merge m:1 sex whoagegrp_`p' using who_`p'_maxvel_3sd.dta, nogen
	drop if subjid==""
	gen who_mindiff_`p'=.
	gen who_maxdiff_`p'=.
	foreach i in 2 3 4 6 {
		replace who_mindiff_`p'=whoinc_`i'_`p' if whoinc_age_`p'==`i'
		replace who_maxdiff_`p'=max_whoinc_`i'_`p' if whoinc_age_`p'==`i'
	}
	replace who_mindiff_`p'=who_mindiff_`p'*d_agedays_`p'/(whoinc_age_`p'*30.4375) if d_agedays_`p'<(whoinc_age_`p'*30.4375)
	replace who_maxdiff_`p'=who_maxdiff_`p'*d_agedays_`p'/(whoinc_age_`p'*30.4375) if d_agedays_`p'>(whoinc_age_`p'*30.4375)
	replace mindiff_`p'=0.5*who_mindiff_`p'-3 if who_mindiff_`p'!=. & d_agedays_`p'<(9*30.4375)
	replace maxdiff_`p'=2*who_maxdiff_`p'+3 if who_maxdiff_`p'!=. & d_agedays_`p'<(9*30.4375)
	replace mindiff_`p'=0.5*who_mindiff_`p'-3 if mindiff_`p'==. & who_mindiff_`p'!=.
	replace maxdiff_`p'=2*who_maxdiff_`p'+3 if maxdiff_`p'==. & who_maxdiff_`p'!=. 
	replace mindiff_`p'=-3 if mindiff_`p'==.
	foreach m in min max{
		sort subjid_`p' agedays
		by subjid_`p': gen `m'diff_prev_`p'=`m'diff_`p'[_n-1]
	}
	*EWMA
	sort subjid_`p' agedays
	by subjid_`p': gen tot_`p'=_N
	by subjid_`p': gen sn_`p'=_n
	summ tot_`p' if exc_`p'==0
	local m=r(max)
	forvalues n=1/`m' {
		by subjid_`p': gen temp_tbc`p'sd_`n'=tbc`p'sd[`n'] 
		by subjid_`p': gen temp_age_`p'_`n'=agedays[`n'] 
		gen double temp_ewt_`p'_`n'=abs(`a' + agedays-temp_age_`p'_`n')^`e' if sn_`p'!=`n'
		gen double temp_ewtsd_`p'_`n'=temp_tbc`p'sd_`n'*temp_ewt_`p'_`n' 
		}
	egen double sumewt_`p'=rowtotal(temp_ewt_`p'_*) if exc_`p'==0
	egen double sumewtsd_`p'=rowtotal(temp_ewtsd_`p'_*) if exc_`p'==0
	gen ewma_`p'=sumewtsd_`p'/sumewt_`p' if subjid_`p'!=""
	gen dewma_`p'=tbc`p'sd-ewma_`p'
	gen d_prev_`p'=.
	sort subjid_`p' agedays
	by subjid_`p': replace d_prev_`p'=`p'-`p'[_n-1] if sn_`p'!=1
	gen d_`p'=.
	by subjid_`p': replace d_`p'=d_prev_`p'[_n+1] if sn_`p'!=tot_`p'	
	gen bef_g_aftm1_`p'=.
	gen aft_g_befp1_`p'=.
	*Before & After
	foreach t in bef aft {
		replace ba="`t'"
		gen double sumewt_`p'_`t'=sumewt_`p'
		gen double sumewtsd_`p'_`t'=sumewtsd_`p'
		gen ewma_`p'_`t'=ewma_`p'
		gen dewma_`p'_`t'=dewma_`p'
		forvalues d=1/`m' {
			if ba=="bef" local s=`d'+1
			if ba=="aft" local s=`d'-1
			replace sumewt_`p'_`t'=sumewt_`p'-temp_ewt_`p'_`d' if temp_ewt_`p'_`d'!=. & sn_`p'==`s'
			replace sumewtsd_`p'_`t'=sumewtsd_`p'-temp_ewtsd_`p'_`d' if temp_ewtsd_`p'_`d'!=. & sn_`p'==`s'
			replace ewma_`p'_`t'=sumewtsd_`p'_`t'/sumewt_`p'_`t' if subjid_`p'!="" & sn_`p'==`s'
			replace dewma_`p'_`t'=tbc`p'sd-ewma_`p'_`t' if subjid_`p'!="" & sn_`p'==`s'
		}
	}	
	sort subjid_`p' agedays
	*Identify which has bigger dewma if other value of measurement pair is not taken into account
	gen pair=1 if (d_prev_`p'<mindiff_prev_`p' | d_`p'<mindiff_`p' | d_prev_`p'>maxdiff_prev_`p'  | d_`p'>maxdiff_`p') & exc_`p'==0
	by subjid_`p': replace bef_g_aftm1_`p'=1 if abs(dewma_`p'_bef)>abs(dewma_`p'_aft[_n-1]) & sn_`p'!=1 & pair==1 & pair[_n-1]==1
	by subjid_`p': replace aft_g_befp1_`p'=1 if abs(dewma_`p'_aft)>abs(dewma_`p'_bef[_n+1]) & sn_`p'!=tot_`p' & pair==1 & pair[_n+1]==1
	
	gen temp_exc_`p'=151 if d_prev_`p'<mindiff_prev_`p' & bef_g_aftm1_`p'==1 & exc_`p'==0 & mindiff_prev_`p'!=.
	gen temp_diff=abs(dewma_`p'_bef) if temp_exc_`p'==151
	
	replace temp_exc_`p'=152 if d_`p'<mindiff_`p' & aft_g_befp1_`p'==1 & exc_`p'==0 & mindiff_`p'!=.
	replace temp_diff=abs(dewma_`p'_aft) if temp_exc_`p'==152
	
	replace temp_exc_`p'=161 if d_prev_`p'>maxdiff_prev_`p' & bef_g_aftm1_`p'==1 & exc_`p'==0 & mindiff_prev_`p'!=.
	replace temp_diff=abs(dewma_`p'_bef) if temp_exc_`p'==161
	
	replace temp_exc_`p'=162 if d_`p'>maxdiff_`p' & aft_g_befp1_`p'==1 & exc_`p'==0 & mindiff_`p'!=.
	replace temp_diff=abs(dewma_`p'_aft) if temp_exc_`p'==162
	
	*Identify more extreme value if only 2 measurements (dewma will be the same)
	gen abstbc`p'sd=abs(tbc`p'sd)
	replace temp_exc_`p'=153 if d_prev_`p'<mindiff_prev_`p' & abstbc`p'sd>abstbc`p'sd[_n-1] & tot_ht==2 & mindiff_prev_`p'!=.	
	replace temp_exc_`p'=154 if d_`p'<mindiff_`p' & abstbc`p'sd>abstbc`p'sd[_n+1] & tot_ht==2 & mindiff_`p'!=.
	replace temp_exc_`p'=163 if d_prev_`p'>maxdiff_prev_`p' & abstbc`p'sd>abstbc`p'sd[_n-1] & tot_ht==2 & maxdiff_prev_`p'!=.
	replace temp_exc_`p'=164 if d_`p'>maxdiff_`p' & abstbc`p'sd>abstbc`p'sd[_n+1] & tot_ht==2 & maxdiff_`p'!=.
	gen temp_bin_exc_`p'=temp_exc_`p'!=.
	gsort subjid_`p' -temp_bin_exc_`p' -temp_diff
	by subjid_`p': gen dn=_n if temp_bin_exc_`p'==1
	replace exc_`p'=int(temp_exc_`p'/10) if dn==1
	replace subjid_`p'="" if exc_`p'>0
	bysort subjid: egen tot_temp_exc_`p'=total(temp_bin_exc_`p')
	egen max_tot_temp_exc_`p'=max(tot_temp_exc_`p')
	local i=max_tot_temp_exc_`p'[1]
	drop ba-max_tot_temp_exc_`p'
}


lab def exc 15 "Min diff" 16 "Max diff", add
save temp_abs_`length', replace
}

************************************************************************************************
* #16: 1 or 2 measurements
************************************************************************************************


foreach length in short medium long {

use temp_abs_`length', clear

foreach p in wt ht {
	sort subjid_`p' agedays
	by subjid_`p': gen vistot_`p'=_N
	by subjid_`p': gen sn_`p'=_n
	gen tbc`p'sd_other=tbc`p'sd[_n+1] if vistot_`p'==2 & sn_`p'==1
	replace tbc`p'sd_other=tbc`p'sd[_n-1] if vistot_`p'==2 & sn_`p'==2
	gen agedays_`p'_other=agedays[_n+1] if vistot_`p'==2 & sn_`p'==1
	replace agedays_`p'_other=agedays[_n-1] if vistot_`p'==2 & sn_`p'==2
	gen absd_tbc`p'sd=abs(tbc`p'sd-tbc`p'sd_other) if vistot_`p'==2
	gen absd_agedays_`p'=abs(agedays-agedays_`p'_other) if vistot_`p'==2
	bysort subjid_`p': egen median_tbc`p'sd_i=median(tbc`p'sd) if exc_`p'==0
	bysort subjid: egen median_tbc`p'sd=min(median_tbc`p'sd_i)
	gen comp_tbc`p'sd_i=tbc`p'sd if exc_`p'==0
	bysort subjid agedays: egen comp_tbc`p'sd=min(comp_tbc`p'sd_i)	
	gen abs_tbc`p'sd=abs(tbc`p'sd) if vistot_`p'==2
}


foreach p in wt ht {
	if `p'==wt local o ht
	if `p'==ht local o wt
	*If only 2 values for a parameter
	sort subjid_`p' agedays
	gen abs_d_`p'_`o'=abs(tbc`p'sd-tbc`o'sd)
	replace abs_d_`p'_`o'=abs(tbc`p'sd-median_tbc`o'sd) if tbc`o'sd==.
	sort subjid_`p' abs_d_`p'_`o' abs_tbc`p'sd 
	by subjid_`p': gen dn_`p'=_n if exc_`p'==0 & vistot_`p'==2
	replace exc_`p'=17 if vistot_`p'==2 & absd_agedays_`p'>365.25 & absd_tbc`p'sd>3 & dn_`p'==2
	replace exc_`p'=18 if vistot_`p'==2 & absd_agedays_`p'<365.25 & absd_tbc`p'sd>2 & dn_`p'==2
	drop vistot_`p'
	bysort subjid_`p': gen vistot_`p'=_N if exc_`p'==0
	*If only 1 value for a parameter
	replace exc_`p'=19 if vistot_`p'==1 & abs(tbc`p'sd)>3 & abs(tbc`p'sd-tbc`o'sd)>5 & tbc`o'sd!=. & exc_`p'==0
	replace exc_`p'=19 if vistot_`p'==1 & abs(tbc`p'sd)>3 & tbc`o'sd==. & abs(tbc`p'sd-median_tbc`o'sd)>5 & median_tbc`o'sd!=. & exc_`p'==0
	replace exc_`p'=19 if vistot_`p'==1 & abs(tbc`p'sd)>5 & tbc`o'sd==. & median_tbc`o'sd==. & exc_`p'==0
}

drop sn_wt-vistot_ht
lab def exc 17 "2 meas >1 year" 18 "2 meas <1 year" 19 "1 meas", add
save temp_1and2_`length', replace

}

************************************************************************************************
* #17: Error load
************************************************************************************************


foreach length in short medium long {

use temp_1and2_`length', replace

foreach p in wt ht {
	*2 measurements
	if `p'==wt local o ht
	if `p'==ht local o wt
	recode exc_`p' 0/1=. 3=. 7=. 20=. 21=. nonm=1, gen(exc_bin_`p')
	recode exc_`p' 0=1 nonm=., gen(inc_bin_`p')
	bysort subjid: egen tot_exc_`p'=total(exc_bin_`p')
	bysort subjid: egen tot_inc_`p'=total(inc_bin_`p')
	*1 measurement
	replace exc_`p'=20 if tot_exc_`p'>=2 & tot_exc_`p'/tot_inc_`p'>0.5 & exc_`p'==0
	replace exc_`o'=21 if tot_exc_`p'>=2 & tot_exc_`p'/tot_inc_`p'>1 & exc_`o'==0
	replace subjid_`p'="" if exc_`p'>0
	replace subjid_`o'="" if exc_`o'>0
}

drop *_bin_* tot_exc_* tot_inc_*

lab def exc 20 "Error load" 21 "Error load other param", add
save temp_load_`length', replace

}



************************************************************************************************
* Combine split datasets
************************************************************************************************

use temp_load_long, clear
append using temp_load_medium
append using temp_load_short
save temp_epros_exclude, replace
save epros_04082015, replace

*Clean up
use temp_epros_exclude, clear
keep subjid-obsid sex exc_wt exc_ht switch replaced* tbcwtsd tbchtsd cdcwtsd cdchtsd cdcwtz cdchtz rcwtsd rchtsd
rename replaced* replace*
drop sex wt ht
rename sexorigcode sex
saveold epros_exclude_04082015, replace

use epros_exclude_04082015, clear
rename exc_*t exc_*t_l
rename switch switch_l
rename replace* replace*_l
rename *sd *sd_l
save epros_exclude_04082015_l, replace










