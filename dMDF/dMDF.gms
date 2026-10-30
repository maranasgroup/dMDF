$INLINECOM   /* */
$onlisting

OPTIONS 
       decimals = 8
       solprint = on
       reslim = 5400
       iterlim = 1000
       domlim = 1000
       limcol = 1000
       limrow = 1000
       optca = 0.0
       optcr = 1E-9
       nlp = ipopt
;


$setglobal nscale 1E3
$setglobal nscaleback 1E-3

SETS
        i 	Set of metabolites 
$include input_files/mets.txt 

        j 	Set of reactions 
$include input_files/rxns.txt

	t	Ordinance of Time Step
	/ 1 * 16 /

	pre_ALDD2x(t)	Time Points Before Aldehyde Dehydrogenase is Active
		/ 1 * 3 /

	pre_PFK(t)	Time Points Before ATP-PFK is Active
		/ 1 * 3 /


	bt(t)	Timepoints with Flux and Biomass Measurements
		/ 1 * 15 /

	ec(i)	Extracellular Metabolites
/
	'ac_e'
	'co2_c'
	'h_c'
	'lac__L_e'
	'mal__L_e'
	'for_e'
	'pyr_e'
	'etoh_c'
	'cell3_e'
	'cell4_e'
	'cell5_e'
	'cell6_e'
	'cellb_e'
	'glc__D_e'
	'h_e'
/
	aux(j)	Auxiliary Reactions

/
	'ND-NAD-FWD'
	'ND-NADH-FWD'
	'ND-FDXR-FWD'
	'ND-FDXO-FWD'
	'ND-NADPH-FWD'
	'ND-NADP-FWD'
	'ND-ADP-ATP-FWD'
	'ND-ATP-ADP-FWD'
	'ND-AMP-ADP-FWD'
	'ND-ADP-AMP-FWD'
	'ND-ATP-AMP-FWD'
	'ND-GTP-FWD'
	'ND-GDP-FWD'
	'ND-PPI-FWD'
	'ND-PI-FWD'
	'ND-PROTON'
	'PH-STABLE'
	'ND-CO2-SEQ'
	'ND-COA-FWD'
	'ND-ACCOA-FWD'
/
	rxn_off(j,t)	
$include input_files/rxn_off.txt
;

PARAMETERS

  S(i,j) 	contains the Stoichiometric matrix of the metabolic model
$include input_files/stoich.txt

  dGS(j) 	Standard Gibbs Free Energy of Reaction

  emf(j)

  Bio_LB(t) 	Biomass Concentration

  Bio_UB(t)	

  tp(t)		Time Step Value

  B_max

  B_max_T(bt)
;

VARIABLES
	        
	c(i,t)      Concentration

	v(j,t)      Flux

	Bio(t)

        z 	    Objective

	DF(j,t)

	B

;

** Set Values At Each Time Point

tp('1') = 8;
tp('2') = 14;
tp('3') = 18;
tp('4') = 24;
tp('5') = 34;
tp('6') = 38;
tp('7') = 43;
tp('8') = 48;
tp('9') = 56;
tp('10') = 62;
tp('11') = 67;
tp('12') = 72;
tp('13') = 82;
tp('14') = 91;
tp('15') = 106;
tp('16') = 115;


** Fractions of Glucose Oligomers Uptaken Per mmol of Cellulose
$setGlobal glc_up_frac 5.20742E-10
$setGlobal cellb_up_frac 2.61383E-05
$setGlobal cell3_up_frac 1.52377E-02
$setGlobal cell4_up_frac 1.61802E-01
$setGlobal cell5_up_frac 6.07311E-02
$setGlobal cell6_up_frac 5.61579E-04

** Concentration Bounds
c.lo(i,t) = .001;
c.up(i,t) = 20;

c.lo(ec,t) = 1E-5;
c.up(ec,t) = 500;

c.up(i,'1') = 20;

c.lo('h_c',t) = 1E-7;
c.up('h_c',t) = 1E-1;

c.fx('h_e',t) = 1E-4;

v.lo(j,t) = 0;
v.up(j,t) = 150;


** Standard Gibbs Free Energies
dGS('ABGPT-FWD') = -17.7;
dGS('CEPA-FWD') = 3.8;
dGS('CEPA2-FWD') = 2.8;
dGS('CEPA3-FWD') = 2.7;
dGS('CEPA4-FWD') = 3.7;
dGS('CEPA5-FWD') = 3.7;
dGS('PGMT-FWD') = -7.1;
dGS('PGI-FWD') = 1.67;
dGS('F6PPT-FWD') = -5.1;
dGS('PFK-FWD') = -14.2;
dGS('FBA-FWD') = 23.9;
dGS('TPI-FWD') = 7.56;
dGS('GAPD-FWD') = 6.3;
dGS('PGK-FWD') = -18.3;
dGS('PGK2-FWD') = -18.54;
dGS('PGM-FWD') = 4.4;
dGS('ENO-FWD') = 1.8;
dGS('PEPCK_re-FWD') = -10.6;
dGS('MDH-FWD') = -30.65;
dGS('ME2-FWD') = 11.35;
dGS('PPDK-FWD') = -19.42;
dGS('LDH_L-FWD') = -27.55;
dGS('PFL-FWD') = -21.02;
dGS('POR-FWD') = -20.9;
dGS('PTAr-FWD') = 9.85;
dGS('ACKr-FWD') = -13.4;
dGS('ALDD2x-FWD') = -46.07;
dGS('ACALD-FWD') = 15.77;
dGS('ALCD2x-FWD') = -21.73;
dGS('PTAr-REV') = -9.85;
dGS('ACKr-REV') = 13.4;
dGS('ACALD-REV') = -15.77;
dGS('ALCD2x-REV') = 21.73;


dGS('EXCH-MAL-FWD') = 0;
dGS('EXCH-MAL-REV') = 0;
dGS('EXCH-PYR-FWD') = 0;
dGS('EXCH-PYR-REV') = 0;
dGS('EXCH-FOR-FWD') = 0;
dGS('EXCH-FOR-REV') = 0;
dGS('EXCH-AC-FWD') = 0;
dGS('EXCH-AC-REV') = 0;
dGS('EXCH-LAC-FWD') = 0;

dGS('EXCH-GLC-FWD') = -30;
dGS('EXCH-CELLB-FWD') = -30;
dGS('EXCH-CELL3-FWD') = -30;
dGS('EXCH-CELL4-FWD') = -30;
dGS('EXCH-CEll5-FWD') = -30;
dGS('EXCH-CELL6-FWD') = -30;

dGS('EXCH-MAL-FWD-ATP') = -30;
dGS('EXCH-PYR-FWD-ATP') = -30;
dGS('EXCH-FOR-FWD-ATP') = -30;
dGS('EXCH-AC-FWD-ATP') = -30;
dGS('EXCH-LAC-FWD-ATP') = -30;

dGS('EXCH-MAL-FWD-H') = 0;
dGS('EXCH-PYR-FWD-H') = 0;
dGS('EXCH-FOR-FWD-H') = 0;
dGS('EXCH-AC-FWD-H') = 0;
dGS('EXCH-LAC-FWD-H') = 0;

emf(j) = 0;

emf('EXCH-MAL-FWD') = -12.26;
emf('EXCH-PYR-FWD') = -12.26;
emf('EXCH-FOR-FWD') = -12.26;
emf('EXCH-AC-FWD') = -12.26;
emf('EXCH-LAC-FWD') = -12.26;

emf('EXCH-MAL-FWD-ATP') = -12.26;
emf('EXCH-PYR-FWD-ATP') = -12.26;
emf('EXCH-FOR-FWD-ATP') = -12.26;
emf('EXCH-AC-FWD-ATP') = -12.26;
emf('EXCH-LAC-FWD-ATP') = -12.26;

emf('EXCH-MAL-REV') = 12.26;
emf('EXCH-PYR-REV') = 12.26;
emf('EXCH-FOR-REV') = 12.26;
emf('EXCH-AC-REV') = 12.26;


B_max = 1000000;

v.fx(rxn_off) = 0;

Bio_LB('1') =0.4462947915;
Bio_LB('2') =0.717938881875;
Bio_LB('3') =1.3489013532;
Bio_LB('4') =3.454798470225;
Bio_LB('5') =5.355164333325;
Bio_LB('6') =5.53699983615;
Bio_LB('7') =5.145095591475;
Bio_LB('8') =4.56782668965;
Bio_LB('9') =4.0240007298;
Bio_LB('10') =3.72166979985;
Bio_LB('11') =3.5747897394;
Bio_LB('12') =3.2364950379;
Bio_LB('13') =3.0837330462;
Bio_LB('14') =2.969959818;
Bio_LB('15') =2.73115367475;
Bio_LB('16') =2.728464792525;


Bio_UB('1') =0.4462947915;
Bio_UB('2') =0.717938881875;
Bio_UB('3') =1.3489013532;
Bio_UB('4') =3.454798470225;
Bio_UB('5') =5.355164333325;
Bio_UB('6') =5.53699983615;
Bio_UB('7') =5.145095591475;
Bio_UB('8') =4.56782668965;
Bio_UB('9') =4.0240007298;
Bio_UB('10') =3.72166979985;
Bio_UB('11') =3.5747897394;
Bio_UB('12') =3.2364950379;
Bio_UB('13') =3.0837330462;
Bio_UB('14') =2.969959818;
Bio_UB('15') =2.73115367475;
Bio_UB('16') =2.728464792525;

$include input_files/conc_bounds.txt

$include input_files/initial_opt.txt

Equations

	Obj		/* Maximization of MDF */

	MDF(j,bt)	/* MDF Lower Bound */

	Dr(j,bt)

$include input_files/TPC_declares.txt	/* Flux-Concentration Constraints */

	FluxLB_AC(pre_ALDD2x)	/* Turn off Aldehyde Dehydrogenase */

	FluxLB_PFK(pre_PFK)	/* Turn off PFK Reaction */

	Bio_L(bt)

	Bio_U(bt)


;

Obj..				z =e= B;


MDF(j,bt)$(not aux(j))..			-1 * (emf(j) + dGS(j) + .008314*333.15 * sum(i, S(i,j) * (log(c(i,bt)) - log(1000)))) =e= DF(j,bt);	

Dr(j,bt)$(not aux(j) and not rxn_off(j,bt))..			DF(j,bt) =g= B;

$include input_files/TPC_eqns.txt

FluxLB_AC(pre_ALDD2x)..		v('ALDD2x-FWD',pre_ALDD2x) =e= 0;

FluxLB_PFK(pre_PFK)..		v('PFK-FWD',pre_PFK) =e= 0;

Bio_L(bt)..				Bio(bt) =g= (Bio_LB(bt) * 0.8);

Bio_U(bt)..				Bio(bt) =l= (Bio_UB(bt) * 1.2);


;

model MaxMin
/

	Obj
	MDF
	Dr
	FluxLB_AC
	FluxLB_PFK
	Bio_L
	Bio_U

$include input_files/TPC_model.txt

/
;

MaxMin.optfile = 1;

Solve MaxMin using nlp maximizing z;

FILE Flux /'results/flux.txt'/;
PUT Flux;


loop(t, 
	Put "Time", system.tab, t.tl/;
	PUT "Reaction", system.tab, "Flux (mmol/gDW)"/;
	loop(j,
		PUT j.tl:0:15, system.tab, v.l(j,t):0:8/;
	);
);


PUTCLOSE;

FILE Conc /'results/conc.txt'/;
PUT Conc;

loop(t, 
	Put "Time", system.tab, t.tl/;
	PUT "Metabolite", system.tab, "Concentration (mmol/L)"/;
	loop(i,
		PUT i.tl, system.tab, (c.l(i,t) * %nscaleback%):0:8/;
	);
);

PUTCLOSE;

FILE Drive /'results/MDF.txt'/;

PUT Drive;

loop(bt, 
	B_max_T(bt) = 10000;
	loop(j$(not aux(j) and not rxn_off(j,bt)),
		if(DF.l(j,bt) < B_max,
			if(v.l(j,bt) > 1E-8,
				B_max = DF.l(j,bt);
			);
		);
		if(DF.l(j,bt) < B_max_T(bt),
			if(v.l(j,bt) > 1E-8,
				B_max_T(bt) = DF.l(j,bt);
			);
		);
	);	
);


PUT "Time", system.tab, "MDF", system.tab, "Rxn"/;

loop(bt,
	loop(j$(not aux(j) and not rxn_off(j,bt)),
		if(DF.l(j,bt) = B_max_T(bt),
			PUT bt.tl, system.tab, B_max_T(bt):0:10, system.tab, j.tl:0:10/;
		);
	);
);

PUT "Min Driving Force", system.tab, B_max:0:10/;

PUTCLOSE;

FILE BioR /'results/bio.txt'/;

PUT BioR;

PUT "Time", system.tab, "Biomass (g/L)"/;

loop(t, 
		PUT t.tl, system.tab, Bio.l(t):0:8/;
);

PUTCLOSE;

FILE thermo_jt / 'results/DF_jt.txt'/;
PUT thermo_jt

loop(j,
	PUT "Reaction",system.tab,j.tl:0:14/;
	PUT "Timepoint",system.tab,"Flux",system.tab,"Driving Force (kJ/mol)"/;	
	loop(t,
		PUT t.tl,system.tab,v.l(j,t):0:10,system.tab,DF.l(j,t):0:10/;
	);
	PUT /;
);

PUTCLOSE; 

FILE thermo_tj / 'results/DF_tj.txt'/;
PUT thermo_tj

loop(t,
	PUT "Timepoint",system.tab,t.tl/;
	PUT "Reaction",system.tab,"Flux",system.tab,"Driving Force (kJ/mol)"/;
	loop(j,
		PUT j.tl:0:14,system.tab,v.l(j,t):0:10,system.tab,DF.l(j,t):0:10/;
	);
	PUT /;
);
 

PUTCLOSE;

FILE Levels / 'input_files/initial_opt.txt'/;

PUT Levels

loop(i,
	loop(t,
		PUT "c.l('",i.tl:0:15,"','",t.tl:0:15,"') = ",(%nscaleback% * c.l(i,t)):0:15,";"/;
	);
);

PUT /;

loop(j,
	loop(t,
		PUT "v.l('",j.tl:0:15,"','",t.tl:0:15,"') = ",v.l(j,t):0:15,";"/;
	);
);

PUT /;


loop(t,
	PUT "Bio.l('",t.tl:0:15,"') = ",Bio.l(t):0:15,";"/;
);


PUT /;

loop(j,
	loop(t,
		PUT "DF.l('",j.tl:0:15,"','",t.tl:0:15,"') = ",DF.l(j,t):0:15,";"/;
	);
);

PUT /;

PUTCLOSE;












