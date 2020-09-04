clear all
set more off
cap log cl
cd "C:\Users\03638881\Dropbox\merror\covid\data\"
grstyle init
grstyle set plain

**************
**XTFRONTIER**
**************
*TI model
*time FEs
*drop $txs

log using covid05.log, replace

qui {
use covid_weekly, clear

xtset id tt2
noi di "Original Sample"
noi distinct id

*Sample* 
g i=1
by id: egen n=total(i)
keep if n==21 
by id: g t=_n
drop if t==1
drop t
sort id tt2
by id: g t=_n
g t2=t^2

noi di "Full Sample: Countries in since 1/1/20"
noi distinct id

*Covariates*
ihstrans tests_new
by id: replace ihs_tests_new=L.ihs_tests_new if _n>1 
replace pop=log(pop)
replace gdppc=log(gdppc)
replace popdens=popdens/1000
glo xs "cpi pop popdens median_age aged_65_older aged_70_older diabetes gdppc"
glo xs2 " "
foreach x in $xs {
	g `x'2=`x'^2	
	glo xs2 "${xs2} `x'2" 
}
glo txs " "
foreach x in $xs {
	g t`x'=t*`x'	
	glo txs "${txs} t`x'" 
}
glo txs2 " "
foreach x in $xs2 {
	g t`x'=t*`x'	
	glo txs2 "${txs2} t`x'" 
}

*Stringency*
g index1 = 0.5*(L.stringencyindex2+L2.stringencyindex2)
g index2 = 0.5*(L3.stringencyindex2+L4.stringencyindex2)
g index3 = 0.33333*(L2.stringencyindex2+L3.stringencyindex2+L4.stringencyindex2)
g index4 = 0.5*(L5.stringencyindex2+L6.stringencyindex2)
g index5 = 0.33333*(L4.stringencyindex2+L5.stringencyindex2+L6.stringencyindex2)

*Sample: Cases*
g r=cases_new>=5
g rt = t if r==1
by id: egen rrt=min(rt)
drop if t<rrt
drop r rt rrt
sort id tt2
by id: g rt=_n

noi di "Full Sample: Obs after 1st week w/ 5+ new cases"
noi distinct id

*Sample: Deaths*
g r=deaths_new>=10
g drt = t if r==1
by id: egen rrt=min(drt)
g dsample=(t>=rrt)
drop r drt rrt
sort id tt2
g rtd=1 if dsample==1 & L.dsample==0
by id: replace rtd=L.rtd + 1 if dsample==1 & rtd==. 

noi di "Full Sample: Obs after 1st week w/ 10+ deaths"
noi distinct id if dsample==1

*Outcomes*
ihstrans cases_new
ihstrans deaths_new
glo ys "ihs_deaths_new ihs_cases_new"

*Scaling*
replace median_age2=median_age2/10

*Final Sample*
keep if cpi!=. & pop!=. & popdens!=. & median_age!=. & aged_65_older!=. & aged_70_older!=. & gdppc!=.
noi su $ys $xs index*, sep(100)

egen reg=group(region)
lab var cpi "Corruption Index"
lab var pop "Log(Population)"
lab var popdens "Population Density"
lab var median_age "Median Age"
lab var aged_65_older "Aged 65+, Percent"
lab var aged_70_older "Aged 70+, Percent"
lab var diabetes "Diabetes, Percent"
lab var gdppc "Log(Per Capita GDP)"

lab var cpi2 "(Corruption Index)$^2$"
lab var pop2 "(Log(Population))$^2$"
lab var popdens2 "(Population Density)$^2$"
lab var median_age2 "(Median Age)$^2$"
lab var aged_65_older2 "(Aged 65+, Percent)$^2$"
lab var aged_70_older2 "(Aged 70+, Percent)$^2$"
lab var diabetes2 "(Diabetes, Percent)$^2$"
lab var gdppc2 "(Log(Per Capita GDP))$^2$"

lab var tcpi "(Corruption Index)x(Linear Time Trend)"
lab var tpop "Log(Population)x(Linear Time Trend)"
lab var tpopdens "(Population Density)x(Linear Time Trend)"
lab var tmedian_age "(Median Age)x(Linear Time Trend)"
lab var taged_65_older "(Aged 65+, Percent)x(Linear Time Trend)"
lab var taged_70_older "(Aged 70+, Percent)x(Linear Time Trend)"
lab var tdiabetes "Diabetes, Percentx(Linear Time Trend)"
lab var tgdppc "Log(Per Capita GDP)x(Linear Time Trend)"

lab var tcpi2 "(Corruption Index)$^2$x(Linear Time Trend)"
lab var tpop2 "(Log(Population))$^2$x(Linear Time Trend)"
lab var tpopdens2 "(Population Density)$^2$x(Linear Time Trend)"
lab var tmedian_age2 "(Median Age)$^2$x(Linear Time Trend)"
lab var taged_65_older2 "(Aged 65+, Percent)$^2$x(Linear Time Trend)"
lab var taged_70_older2 "(Aged 70+, Percent)$^2$x(Linear Time Trend)"
lab var tdiabetes2 "(Diabetes, Percent)$^2$x(Linear Time Trend)"
lab var tgdppc2 "(Log(Per Capita GDP))$^2$x(Linear Time Trend)"

lab var index1 "0.5*(Index$_{t-1}$+Index$_{t-2}$"
lab var index2 "0.5*(Index$_{t-3}$+Index$_{t-4}$"
lab var index3 "0.33*(Index$_{t-2}$+Index$_{t-3}$+Index$_{t-4}$"
lab var index4 "0.5*(Index$_{t-5}$+Index$_{t-6}$"
lab var index5 "0.33*(Index$_{t-4}$+Index$_{t-5}$+Index$_{t-6}$"

*Final Covariates
glo xs "cpi pop popdens median_age aged_65_older aged_70_older diabetes gdppc"
glo xs "pop popdens median_age aged_65_older gdppc cpi "
glo xs2 " "
foreach x in $xs {	
	glo xs2 "${xs2} `x'2" 
}
glo txs " "
foreach x in $xs {
	glo txs "${txs} t`x'" 
}
glo txs2 " "
foreach x in $xs2 {
	glo txs2 "${txs2} t`x'" 
}

glo xsave "${xs} ${xs2}"
glo xs "${xs} ${xs2} i.reg i.rt"
constraint 1 [mu]_cons=0

**************************
**Estimation: xtfrontier**
**************************

loc y "ihs_cases_new"
noi di "**** `y' ****"

forval i=1/3 {

	xtfrontier `y' $xs index`i', ti constraint(1)
	eststo m`i'
	predict u, u
	g z=`y' if e(sample)
	g y=z+u
	g yy=0.5*(exp(y) - exp(-y))
	g zz=0.5*(exp(`y') - exp(-`y'))
	replace zz=. if yy==.
	egen s=total(zz)
	egen ss=total(yy)
	su s
	eststo m`i', add(Sum r(mean))
	su ss
	eststo m`i', add(SFASum r(mean))
	if (`i'==2) {
	    preserve
		collapse (sum) zz yy (mean) pop, by(iso country)
		replace pop=exp(pop)
		g irate1 = 100*zz/pop
		g irate2 = 100*yy/pop
		sort country
		noi l iso country irate1 irate2, sep(100)
		restore
		lab var zz "New Cases, Observed"
		lab var yy "New Cases, Predicted"
		noi tsline zz if iso=="USA", ms(i) lc(gs0) lp(solid) || tsline yy if iso=="USA", ms(i) lc(gs6) lp(dash) ///
			subtitle("United States") saving(grxt0, replace) ylab(,labs(small) angle(0)) xlab(,labs(small)) xtitle(" ") ///
			legend(col(1) size(small) rowg(*.5) keyg(*.5) symysize(*.5) symxsize(*.8))
		graph export "C:\Users\03638881\Dropbox\merror\covid\text\grxt20.png", replace
		preserve
		keep id iso country zz yy region tt2
		save xtcases, replace
		collapse (sum) zz yy, by(region tt2)
		egen r=group(region)
		lab var zz "Aggregate New Cases, Observed"
		lab var yy "Aggregate New Cases, Predicted"
		forval s=1/7 {
			if (`s' == 1) loc note "East Asia & Pacific"
			if (`s' == 2) loc note "Europe & Central Asia"
			if (`s' == 3) loc note "Latin America & Caribbean"
			if (`s' == 4) loc note "Middle East & North Africa"
			if (`s' == 5) loc note "North America"
			if (`s' == 6) loc note "South Asia"
			if (`s' == 7) loc note "Sub-Saharan Africa"
			noi tsline zz if r==`s', ms(i) lc(gs0) lp(solid) || tsline yy if r==`s', ms(i) lc(gs6) lp(dash) ///
				subtitle("`note'") saving(grxt`s', replace) ylab(,labs(small) angle(0)) xlab(,labs(small)) xtitle(" ") ///
				legend(col(1) size(small) rowg(*.5) keyg(*.5) symysize(*.5) symxsize(*.8))
			graph export "C:\Users\03638881\Dropbox\merror\covid\text\grxt2`s'.png", replace
		}
		restore
	}
	drop s ss z zz y yy u

}

noi estout m1 m2 m3, k(*) cells(b(star fmt(3)) se(nopar fmt(3))) stats(N Sum SFASum, fmt(%9.0g 0 0)) style(fixed) starlevels(a 0.10 b 0.05 c 0.01) legend

loc y "ihs_deaths_new"
noi di "**** `y' ****"

if ("`y'" == "ihs_deaths_new") loc sample "if dsample==1"

forval i=3/5 {
	
	xtfrontier `y' $xs index`i' `sample', ti constraint(1)
	eststo d`i'
	predict u, u
	g z=`y' if e(sample)
	g y=z+u
	g yy=0.5*(exp(y) - exp(-y))
	g zz=0.5*(exp(`y') - exp(-`y'))
	replace zz=. if yy==.
	egen s=total(zz)
	egen ss=total(yy)
	su s
	eststo d`i', add(Sum r(mean))
	su ss
	eststo d`i', add(SFASum r(mean))
	if (`i'==4) {
		lab var zz "Weekly Deaths, Observed"
		lab var yy "Weekly Deaths, Predicted"
		noi tsline zz if iso=="USA", ms(i) lc(gs0) lp(solid) || tsline yy if iso=="USA", ms(i) lc(gs6) lp(dash) ///
			subtitle("United States") saving(grdxt0, replace) ylab(,labs(small) angle(0)) xlab(,labs(small)) xtitle(" ") ///
			legend(col(1) size(small) rowg(*.5) keyg(*.5) symysize(*.5) symxsize(*.8))
		graph export "C:\Users\03638881\Dropbox\merror\covid\text\grdxt20.png", replace
		preserve
		keep id iso country zz yy region tt2
		save xtdeaths, replace
		collapse (sum) zz yy, by(region tt2)
		egen r=group(region)
		lab var zz "Aggregate Weekly Deaths, Observed"
		lab var yy "Aggregate Weekly Deaths, Predicted"
		forval s=1/7 {
			if (`s' == 1) loc note "East Asia & Pacific"
			if (`s' == 2) loc note "Europe & Central Asia"
			if (`s' == 3) loc note "Latin America & Caribbean"
			if (`s' == 4) loc note "Middle East & North Africa"
			if (`s' == 5) loc note "North America"
			if (`s' == 6) loc note "South Asia"
			if (`s' == 7) loc note "Sub-Saharan Africa"
			noi tsline zz if r==`s', ms(i) lc(gs0) lp(solid) || tsline yy if r==`s', ms(i) lc(gs6) lp(dash) ///
				subtitle("`note'") saving(grdxt`s', replace) ylab(,labs(small) angle(0)) xlab(,labs(small)) xtitle(" ") ///
				legend(col(1) size(small) rowg(*.5) keyg(*.5) symysize(*.5) symxsize(*.8))
			graph export "C:\Users\03638881\Dropbox\merror\covid\text\grdxt2`s'.png", replace
		}
		restore
		preserve 
		use xtdeaths, clear
		ren id id2
		merge 1:1 iso tt2 using mortality
		g sample=_m==3
		sort id tt2
		by id: egen s=sum(sample)
		keep if s>=3
		drop _m sample s
		keep if year>=2017
		collapse (sum) dtotal zz yy, by(year week)
		sort year week
		egen tt2=group(year week)
		tsset tt2
		replace dtotal=. if year==2020 & week>18
		g em=dtotal-L52.dtotal if year==2020
		noi tw connected dtotal week if year==2017 & week<20, ms(i) lp(solid) lc(gs10) sort || ///
			connected dtotal week if year==2018 & week<20, ms(i) lp(solid) lc(gs5) sort || ///
			connected dtotal week if year==2019 & week<20, ms(i) lp(solid) lc(gs0) sort || ///
			connected dtotal week if year==2020 & week<20, ms(i) lp(solid) lc(maroon) sort || ///
			connected zz week if year==2020 & week<20, ms(i) lp(shortdash) lc(navy) sort || ///
			connected yy week if year==2020 & week<20, ms(i) lp(longdash_dot) lc(dkgreen) sort || ///
			connected em week if year==2020 & week<20, ms(i) lp(shortdash_dot) lc(dkorange) sort ///
			subtitle("Aggregate Weekly Deaths") saving(grxtem, replace) ylab(,labs(small) angle(0)) xlab(,labs(small)) xtitle("Calendar Week") ///
			legend(lab(1 "All-Cause Mortality, 2017") lab(2 "All-Cause Mortality, 2018") lab(3 "All-Cause Mortality, 2019") lab(4 "All-Cause Mortality, 2020") ///
				lab(5 "Covid-19 Deaths, Observed, 2020") lab(6 "Covid-19 Deaths, SFA, 2020") lab(7 "Covid-19 Deaths, EM, 2020") col(2) colfirst size(small) rowg(*.5) keyg(*.5) symysize(*.5) symxsize(*.8))
		gr export "C:\Users\03638881\Dropbox\merror\covid\text\grxtem2.png", replace
		bys year: g c_dtotal=sum(dtotal)
		bys year: g c_zz=sum(zz)
		bys year: g c_yy=sum(yy)
		bys year: g c_em=sum(em)
		noi tw connected c_dtotal week if year==2017 & week<20, ms(i) lp(solid) lc(gs10) sort || ///
			connected c_dtotal week if year==2018 & week<20, ms(i) lp(solid) lc(gs5) sort || ///
			connected c_dtotal week if year==2019 & week<20, ms(i) lp(solid) lc(gs0) sort || ///
			connected c_dtotal week if year==2020 & week<20, ms(i) lp(solid) lc(maroon) sort || ///
			connected c_zz week if year==2020 & week<20, ms(i) lp(shortdash) lc(navy) sort || ///
			connected c_yy week if year==2020 & week<20, ms(i) lp(longdash_dot) lc(dkgreen) sort || ///
			connected c_em week if year==2020 & week<20, ms(i) lp(shortdash_dot) lc(dkorange) sort ///
			subtitle("Cumulative Weekly Deaths") saving(grxtemc, replace) ylab(,labs(small) angle(0)) xlab(,labs(small)) xtitle("Calendar Week") ///
			legend(lab(1 "All-Cause Mortality, 2017") lab(2 "All-Cause Mortality, 2018") lab(3 "All-Cause Mortality, 2019") lab(4 "All-Cause Mortality, 2020") ///
				lab(5 "Covid-19 Deaths, Observed, 2020") lab(6 "Covid-19 Deaths, Predicted, 2020") lab(7 "Covid-19 Deaths, EM, 2020") col(2) colfirst size(small) rowg(*.5) keyg(*.5) symysize(*.5) symxsize(*.8))
		gr export "C:\Users\03638881\Dropbox\merror\covid\text\grxtemc2.png", replace
		restore
	}
	drop s ss z zz y yy u

}

noi estout d3 d4 d5, k(*) cells(b(star fmt(3)) se(nopar fmt(3))) stats(N Sum SFASum, fmt(%9.0g 0 0)) style(fixed) starlevels(a 0.10 b 0.05 c 0.01) legend

/*
esttab m1 m2 m3 d3 d4 d5 ///
	using "C:\Users\03638881\Dropbox\merror\covid\text\xt2.tex", ///
	se nonotes  style(tex)  b(%12.3f) se(%12.3f)  noobs ///
	starlevels(* 0.10 ** 0.05 *** 0.01) label mlabels("" "" "" "" "" "" "" "" "" "" ) ///
	 keep(${xsave} index*) nonumbers  replace  fragment  ///
	 scalar(N Sum SFASum) prehead({ \begin{tabular}{lcccccc}    ///
	 \hline \hline ///
	 & \multicolumn{3}{c}{Cases} & \multicolumn{3}{c}{Deaths} \\  ///
	 & (1) & (2) & (3) & (4) & (5) & (6) \\ ) ///
	 prefoot( \\ ) ///
	 postfoot(\hline \hline \end{tabular} } \begin{tablenotes}[para,flushleft] \footnotesize{Standard errors in parentheses. ///
	 Time and region fixed effects included in all models. * p $<$.10, ** p$<$ .05, *** p$<$.01.} \end{tablenotes} )
*/


}
log cl
