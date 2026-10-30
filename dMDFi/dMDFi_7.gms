$INLINECOM   /* */
$onlisting

OPTIONS 
       decimals = 8
       solprint = on
       reslim = 5400
       iterlim = 100000
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

	t_fx(t)		What We Fix

	/1*6/

	t_pres(t)	What We Write Down
	/1*7/


	n	Number of Starts
	/ 1 /

	pre_ALDD2x(t)	Time Points Before Aldehyde Dehydrogenase is Active
		/ 1 * 3 /

	pre_PFK(t)	Time Points Before ATP-PFK is Active
		/ 1 * 3 /


	bt(t)	Timepoints with Flux Bounds
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
	aux(j)	Reactions I don't care about for dG purposes

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

c.lo('glc__D_c',t) = 1E-7;
c.lo('cellb_c',t) = 1E-7;
c.lo('cell3_c',t) = 1E-7;
c.lo('cell4_c',t) = 1E-7;
c.lo('cell5_c',t) = 1E-7;
c.lo('cell6_c',t) = 1E-7;

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

$include input_files/res/t6.txt

Equations

	Obj		/* Maximization of MDF */

	MDF(j,bt)	/* MDF Lower Bound */

	Dr(j,bt)

	TPC(i)

	Bio_L(bt)

	Bio_U(bt)


;

Obj..				z =e= B;

MDF(j,'7')$(not aux(j))..			-1 * (emf(j) + dGS(j) + .008314*333.15 * sum(i, S(i,j) * (log(c(i,'7')) - log(1000)))) =e= DF(j,'7');	

Dr(j,'7')$(not aux(j) and not rxn_off(j,'7'))..			DF(j,'7') =g= B;

TPC(i)..	(c(i,'7')) - (c(i,'6')) =e= sum(j, S(i,j) * (v(j,'6') * Bio('6'))) * (tp('7') - tp('6'));

Bio_L('6')..				Bio('6') =g= (Bio_LB('6') * 0.8);

Bio_U('6')..				Bio('6') =l= (Bio_UB('6') * 1.2);

;

model MaxMin
/

	Obj
	MDF
	Dr
	TPC
	Bio_L
	Bio_U
/
;

Solve MaxMin using nlp maximizing z;

FILE prep_t7 /'input_files/res/t7.txt'/;

PUT prep_t7;

loop(t_pres,
	loop(i,
		PUT "c.fx('",i.tl:0:10,"','",t_pres.tl:0:10,"') = ",c.l(i,t_pres):0:10,";"/;
	);
);

PUT /;

loop(t_fx,
	loop(j,
		PUT "v.fx('",j.tl:0:10,"','",t_fx.tl:0:10,"') = ",v.l(j,t_fx):0:10,";"/;
	);
);

PUT /;

loop(t_fx,
	PUT "Bio.fx('",t_fx.tl:0:10,"') = ",Bio.l(t_fx):0:10,";"/;
);
	
PUTCLOSE;

FILE res_t7 /'results/t7.txt'/;

PUT res_t7;

loop(t_pres,
	loop(i,
		PUT "c('",i.tl:0:10,"','",t_pres.tl:0:10,"') = ",c.l(i,t_pres):0:10,";"/;
	);
);

PUT /;

loop(t_fx,
	loop(j,
		PUT "v('",j.tl:0:10,"','",t_fx.tl:0:10,"') = ",v.l(j,t_fx):0:10,";"/;
	);
);

loop(t_pres,
	loop(j,
		PUT "DF('",j.tl:0:10,"','",t_pres.tl:0:10,"') = ",DF.l(j,t_pres):0:10,";"/;
	);
);


PUT /;

loop(t_fx,
	PUT "Bio('",t_fx.tl:0:10,"') = ",Bio.l(t_fx):0:10,";"/;
);

	
PUTCLOSE;





