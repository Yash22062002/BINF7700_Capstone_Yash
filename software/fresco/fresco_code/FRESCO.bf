/* 
FRESCO
Finding Regions of Excess Synonymous Constraint
Rachel Sealfon, 2014
Some functions adapted from HYPHY batch scripts by Sergei Kosokovsky-Pond 
*/


/* defined from codons.bf by Kosokovsky-Pond */
UniversalGeneticCode = {			
	{14,/*AAA*/ 13,/*AAC*/ 14,/*AAG*/  13,/*AAT*/
	7, /*ACA*/ 7, /*ACC*/ 7, /*ACG*/  7, /*ACT*/
	19, /*AGA*/ 5, /*AGC*/19, /*AGG*/  5, /*AGT*/
	2, /*ATA*/ 2, /*ATC*/	3, /*ATG*/  2, /*ATT*/
	12,/*CAA*/ 11,/*CAC*/ 12,/*CAG*/  11,/*CAT*/
	6, /*CCA*/ 6, /*CCC*/ 6, /*CCG*/  6, /*CCT*/
	19,/*CGA*/ 19,/*CGC*/ 19,/*CGG*/  19,/*CGT*/
	1, /*CTA*/ 1, /*CTG*/ 1, /*CTC*/  1, /*CTT*/
	16,/*GAA*/ 15,/*GAC*/ 16,/*GAG*/  15,/*GAT*/
	8, /*GCA*/ 8, /*GCC*/ 8, /*GCG*/  8, /*GCT*/
	20,/*GGA*/ 20,/*GGC*/ 20,/*GGG*/  20,/*GGT*/
	4, /*GTA*/ 4, /*GTC*/ 4, /*GTG*/  4, /*GTT*/
	10,/*TAA*/  9, /*TAC*/10,/*TAG*/   9, /*TAT*/
	5, /*TCA*/ 5, /*TCC*/ 5, /*TCG*/  5, /*TCT*/
	10,/*TGA*/ 17,/*TGC*/ 18,/*TGG*/  17,/*TGT*/		  
	1, /*TTA*/ 0, /*TTC*/ 1, /*TTG*/  0  /*TTT*/ }	    
	};

/*function from codons.bf by Kosokovsky-Pond*/

function BuildCodonFrequencies (nucFreqMatrix){	
	PIStop = 1.0; 		/* denominator */	
	result = {61,1};    /* resulting codon frequencies */	
	hshift = 0;         /* how many stop codons have been counted so far */	
	for (h=0; h<64; h=h+1) /* loop over all possible codons */	{	
		first  = h$16;    /* Decompose a codon into 3 nucleotides.
	 			The index of the first nucleotide (A=0,C=1,G=2,T=3) is found here,
				by doing integer division by 16  */
		second = h%16$4;  /* The index of the second nucleotide.
 				First take the remainder of division by 16, i.e. positions 2 and 3
				and then extract position 2 by integer division by 4*/		
		third  = h%4;     /* The index of the third nucleotide.						
				 Remainder of integer division by 4*/
				/* in the end: h = 16*first + 4*second + third */
		if (UniversalGeneticCode[h]==10) /* stop codon */ {
			hshift = hshift+1; 			
			PIStop = PIStop-nucFreqMatrix[first][0]*nucFreqMatrix[second][1]*nucFreqMatrix[third][2]; /* adjust the denominator */
		}		
		else {			
			result[h-hshift] = nucFreqMatrix[first][0]*nucFreqMatrix[second][1]*nucFreqMatrix[third][2];
 				/* store the frequency for codon h. Notice the substraction of hshift to compensate
				for the absense of stop codons. The first codon affected by it is
 				TAC (h=49), which gets stored in result[48], because TAA (a stop codon) was skipped. */	
		}	
	}	
	return result*(1.0/PIStop);
}

/*Nucleotide model fit*/

REPLACE_TREE_STRUCTURE = 1;
SetDialogPrompt("Please specify a codon data file:");
DataSet myData = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter nucData = CreateFilter (myData, 1);

HarvestFrequencies(overallFrequencies, nucData, 1,1,1);
global kappa = 1;
HKY85RateMatrix = {{*,t,kappa*t,t}{t,*,t,kappa*t}{kappa*t,t,*,t}{t,kappa*t,t,*}};
Model HKY85 = (HKY85RateMatrix, overallFrequencies);


ExecuteAFile(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "queryTree.bf");
Tree myTree = treeString;

LikelihoodFunction nucLF = (nucData,myTree);
Optimize(res,nucLF);

/*Codon model fit*/
DataSetFilter codonFilter = CreateFilter (myData,3, "", "", "TAA,TAG,TGA");
HarvestFrequencies (nuc3by4,codonFilter,3,1,1); /*collect position-specific frequencies*/
estimatedCodonFreqs = BuildCodonFrequencies (nuc3by4);
global R = 1;  /* global dN/dS ratio */

MG94xHKY85RateMatrix = {61,61};

hshift = 0;

for (h=0; h<63; h=h+1){
	if (UniversalGeneticCode[h] == 10) /* stop codon */{
		hshift = hshift + 1;	
	}
	else {
		vshift = hshift;		
		for (v=h+1; v<64; v=v+1) {
			if (UniversalGeneticCode[v] == 10)  /* stop codon */ {
				vshift = vshift + 1;
			}
			else {
				diff = v-h;
				if ((h$4==v$4)||((diff%4==0)&&(h$16==v$16))||(diff%16==0)) /* differ by one subsitution only */ {
					if (h$4==v$4) /* third position */ {
						transition = v%4; /* get targets nucleotides for h->v and v->h substitutions */					
						transition2= h%4;		
						nucPosInCodon = 2; /* change is in the third position */	
					}
					else {					
						if(diff%16==0) /* first position */ {	
							transition = v$16; /* get targets nucleotides for h->v and v->h substitutions */
							transition2= h$16;	
							nucPosInCodon = 0; /* change is in the first position */	
						}					
						else {						
							transition = v%16$4; /* get targets nucleotides for h->v and v->h substitutions */
							transition2= h%16$4;							
							nucPosInCodon = 1;  /* change is in the second position */		
						}					
					}
					if (UniversalGeneticCode[h]==UniversalGeneticCode[v]) /* synonymous */ {
						if (Abs(transition-transition2)%2) /* transversion: difference is not divisible by 2 */	{
							MG94xHKY85RateMatrix[h-hshift][v-vshift] := codonrate*nuc3by4__[transition__][nucPosInCodon__];	
							MG94xHKY85RateMatrix[v-vshift][h-hshift] := codonrate*nuc3by4__[transition2__][nucPosInCodon__];
			  			}
						else /* transition */ {				
							MG94xHKY85RateMatrix[h-hshift][v-vshift] := codonrate*kappa__*nuc3by4__[transition__][nucPosInCodon__];		
							MG94xHKY85RateMatrix[v-vshift][h-hshift] := codonrate*kappa__*nuc3by4__[transition2__][nucPosInCodon__];	
 						}
					}
					else /* non-synonymous */ {
						if (Abs(transition-transition2)%2) /* transversion: difference is not divisible by 2 */	{	
							MG94xHKY85RateMatrix[h-hshift][v-vshift] := R*codonrate*nuc3by4__[transition__][nucPosInCodon__];
							MG94xHKY85RateMatrix[v-vshift][h-hshift] := R*codonrate*nuc3by4__[transition2__][nucPosInCodon__];
			  			}
				
						else /* transition */ {						
							MG94xHKY85RateMatrix[h-hshift][v-vshift] := R*kappa__*codonrate*nuc3by4__[transition__][nucPosInCodon__];
							MG94xHKY85RateMatrix[v-vshift][h-hshift] := R*kappa__*codonrate*nuc3by4__[transition2__][nucPosInCodon__];
			  			}
					}
				}
			}
		}
	}
}
Model MG94xHKY85 = (MG94xHKY85RateMatrix,estimatedCodonFreqs,0);
Tree codonTree = treeString;
LikelihoodFunction lf = (codonFilter, codonTree);
Optimize (myRes, lf);
/*
fprintf (stdout, "\n", lf, "\n");

fprintf (stdout, "\n", kappa, "\n");

fprintf (stdout, "\n", SG, "\n");
fprintf (stdout, "\n", NSG, "\n");
*/

fscanf(stdin, "Number", w);

/*Now, define and optimise a site-by-site matrix.*/

global S = 1; 
global NS = 1;

SitewiseRateMatrix = {61,61};

hshift = 0;

for (h=0; h<63; h=h+1){
	if (UniversalGeneticCode[h] == 10) /* stop codon */{
		hshift = hshift + 1;	
	}
	else {
		vshift = hshift;		
		for (v=h+1; v<64; v=v+1) {
			if (UniversalGeneticCode[v] == 10)  /* stop codon */ {
				vshift = vshift + 1;
			}
			else {
				diff = v-h;
				if ((h$4==v$4)||((diff%4==0)&&(h$16==v$16))||(diff%16==0)) /* differ by one subsitution only */ {
					if (h$4==v$4) /* third position */ {
						transition = v%4; /* get targets nucleotides for h->v and v->h substitutions */					
						transition2= h%4;		
						nucPosInCodon = 2; /* change is in the third position */	
					}
					else {					
						if(diff%16==0) /* first position */ {	
							transition = v$16; /* get targets nucleotides for h->v and v->h substitutions */
							transition2= h$16;	
							nucPosInCodon = 0; /* change is in the first position */	
						}					
						else {						
							transition = v%16$4; /* get targets nucleotides for h->v and v->h substitutions */
							transition2= h%16$4;							
							nucPosInCodon = 1;  /* change is in the second position */		
						}					
					}
					if (UniversalGeneticCode[h]==UniversalGeneticCode[v]) /* synonymous */ {
						if (Abs(transition-transition2)%2) /* transversion: difference is not divisible by 2 */	{
							SitewiseRateMatrix[h-hshift][v-vshift] := S*codonrate*nuc3by4__[transition__][nucPosInCodon__];	
							SitewiseRateMatrix[v-vshift][h-hshift] := S*codonrate*nuc3by4__[transition2__][nucPosInCodon__];
			  			}
						else /* transition */ {				
							SitewiseRateMatrix[h-hshift][v-vshift] := S*codonrate*kappa__*nuc3by4__[transition__][nucPosInCodon__];		
							SitewiseRateMatrix[v-vshift][h-hshift] := S*codonrate*kappa__*nuc3by4__[transition2__][nucPosInCodon__];	
 						}
					}
					else /* non-synonymous */ {
						if (Abs(transition-transition2)%2) /* transversion: difference is not divisible by 2 */	{	
							SitewiseRateMatrix[h-hshift][v-vshift] := NS*R__*codonrate*nuc3by4__[transition__][nucPosInCodon__];
							SitewiseRateMatrix[v-vshift][h-hshift] := NS*R__*codonrate*nuc3by4__[transition2__][nucPosInCodon__];
			  			}
				
						else /* transition */ {						
							SitewiseRateMatrix[h-hshift][v-vshift] := NS*R__*codonrate*kappa__*nuc3by4__[transition__][nucPosInCodon__];
							SitewiseRateMatrix[v-vshift][h-hshift] := NS*R__*codonrate*kappa__*nuc3by4__[transition2__][nucPosInCodon__];
			  			}
					}
				}
			}
		}
	}
}

fprintf(stdout, "#window_start\tsyn_rate_alternate_model\tsyn_rate_null_model\tnonsynonymous_rate_alternate_model\tnonsynonymous_rate_null_model\tlog_likelihood_alternate_model\tlog_likelihood_null_model\tn_local_params_alternate_model\tn_local_params_null_model\tLRT\tp_value\n");
for (siteCount = 0; siteCount < myData.sites/3 - w + 1; siteCount = siteCount+1) {	
	
	filterString = "" + (siteCount*3) + "-" + (siteCount*3+w*3-1); /*pick out individual codons*/
	DataSetFilter siteFilter = CreateFilter (myData,3,filterString,"","TAA,TAG,TGA");
	Model sitewise = (SitewiseRateMatrix,estimatedCodonFreqs,0);
	Tree siteTree = treeString;
	
	ReplicateConstraint ("this1.?.codonrate:=this2.?.codonrate__",siteTree,codonTree);	
	LikelihoodFunction siteLikelihood = (siteFilter, siteTree);
	/*fprintf(siteFilter, "\n");
	fprintf(stdout, "Site #:", siteCount, "\n");*/
	S = 1;
	Optimize(site_res, siteLikelihood);
	lnLA = site_res[1][0];
	dfA = site_res[1][1];
	SA = S;
	NA = NS;
	S := 1;
	Optimize(site_res, siteLikelihood);
	/*fprintf(stdout, siteLikelihood);*/
	lnL0 = site_res[1][0];
	df0 = site_res[1][1];
	LRT = -2*(lnL0-lnLA);
	Pvalue=1-CChi2(LRT,1);
	/*fprintf(stdout,"\n\nThe statistic ",LRT," has P-value ", Pvalue,"\n\n");*/
	fprintf(stdout,siteCount, "\t", SA, "\t", S, "\t", NA, "\t", NS, "\t", lnLA, "\t", lnL0, "\t", dfA, "\t", df0, "\t", LRT,"\t", Pvalue,"\n");
}
