/*
 * data.cpp
 *
 *      Author: Mironov
 */


#include "track_util.h"


const char* version="2.40";




Chromosome *chrom_list;       // list of chromosomes
Chromosome *curChrom=chrom_list;
int  binSize=100;   // frame size fo profile
bool  NAFlag=0;


long long GenomeLength=0;      // TOTAL LENGTH OF THE GENOME
int n_chrom;
char trackName[4096];       // current track name
char *chromFile=0;
char *confFile=0;		// confounder file name
char *cfgFile=0;		// config file name
char *profPath=strdup("./");
char *trackPath=strdup("./");
char *resPath=strdup("./");
char *statFileName=(char*)"./statistics";
char *paramsFileName=(char*)"./params";
char *inputProfiles=0;
AliasTable *aliases=0;
char* aliasFile=0;
char* Rscript=0;
int  plotH=3;
int  plotW=7;



//char *outTrackFile=0; // Filename for write out track
//char *idSuff=(char*)"";


bool  syntax=1;				// Strong syntax control


//int   inpThreshold=0;			// Testing of binarized input data, % of max
int   complFg=IGNORE_STRAND;
int   profileLength;			// size of the profile array


//char *pcorProfile=0;    		// partial correlation profile file name
int  binBufSize=30000000;




int 	kernelType=KERN_NORM;
char* 	customKern=0;
double 	noiseLevel=0;
int 	wSize=100000;        	// size of widow (nucleotides)
int 	wStep=0;             	// window step   (nucleotides)
int 	flankSize=0;
double 	kernelSigma=1000.;    	// kernel width (nucleotides)
double 	kernelShift=0;      	// Kernel mean (for Gauss) or Kernel start for exponent
int 	intervFg0;
double 	scaleFactor=0.2;


//==================== Output
int   	writeDistr=DISTR_DETAIL;
bool  	writeDistCorr=1;		    // write BroadPeak
int   	crossWidth=10000;
bool  	outSpectr=0;
bool  	outChrom=0;
int   	outRes=XML|TAB;
double 	totCorr=FNA, BgTotal=FNA;
int 	RScriptFg=0;
//bool 	writeHTML=false;
//bool  	writePDF=false;


//===================================== Local correlation track
int 	outLC=0;
//int   	lcFlag=CENTER;
int 	LCScale=LOG_SCALE;
double 	L_LC=-20;		// Left treshold to write the Local Correlation track
double 	R_LC=1;			// Right treshold on write the Local Correlation track


//=================================================================
bool 	outPrjBGr=true;


int 	wProfStep=0;          	// window step   (profile scale)
int 	wProfSize=0;          	// size of widow (profile scale)
int 	LFlankProfSize=0;         // size of flank (profile scale)
int 	RFlankProfSize=0;         // size of flank (profile scale)
int 	profWithFlanksLength=0; 	// size of profWindow array (including random flanks)
bool	localSuffle=0;			// use shuffle inside the windoww
double 	kernelProfSigma=1000;     // kernel width ((profile scale)
double 	kernelProfShift=0;
double 	kernelNS=0;			// Correction for non-specifisity
Track 	*track1=0, *track2=0, *projTrack=0;
Kernel 	*kern=0;
double 	maxNA0=95;
double 	maxZero0=95;
double 	maxNA=100;
double 	maxZero=100;
int 	nShuffle=10000;
char 	*trackName1=strdup("");
char 	*trackName2=strdup("");
double 	mannW_Z=0;
double 	mannW_p=1;
double 	smoothZ=3;


Model 	*model;


int 	threshold=0;


FILE 	*logFile=0;
bool 	doAutoCorr=0;


int 	corrScale=10;
double 	prod11=0,prod12=0,prod22=0, eprod1,eprod2;
int 	nprod=0;
XYCorrelation XYfgCorrelation;		    // array for correlation picture
XYCorrelation XYbgcorrelation;			// array for correlation picture
Fourier LCorrelation;


int 	bpType=BP_SIGNAL;
int 	cage=0;
bool 	clearProfile=false;
int 	scoreType=AV_SCORE;
FileListEntry files[256];
int   	nfiles=0;
bool LCExists=false;


//double 	BgAvCorr=0;
//double 	FgAvCorr=0;
int  	pgLevel=2;
float total=0;						// total count over the track
