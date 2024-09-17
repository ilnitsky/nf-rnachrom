#include "track_util.h"
//==============================================================================
// for debugging :
// set debugFg=DEBUG_LOG|DEBUG_PRINT
// set debS string for module identification//
// Use:  deb(n);  // print debug information as number
//       deb(format,...)    // print debug information as printf
//       deb(n,format,...)  // print debug information as a number and printf
// example:
// debS="fun1";
// deb(1);
// ....
// deb(2,"%i %f", n, d);
// ....
// deb("OK");


const char * progName="StereoGene";
const int progType=SG;




void printMiniHelp(){
	printf("\n");
	printf("The %s program compares pairs of tracks and calculates kernel correlations\n",progName);
	printf("===========  version %s ========\n",version);
	printf("Usage:\n");
	printf("$ ./%s [-parameters] trackFile_1 trackFile_2 ... trackFile_n\n",progName);
	printf("\n");
	printf("Say %s -h for more information\n",progName);
	printf("\n");
	exit(0);
}


void printProgDescr(){
	printf("\n");
	printf("The StereoGene program compares pairs of tracks and calculates kernel correlations\n");
	printf("Usage:\n");
	printf("$ ./StereoGene [-parameters] trackFile_1 trackFile_2 ... trackFile_n\n");
	printf("\n");
}


//============================================ Tests =========================================
#include <unistd.h>
void tst(){
	debugFg=DEBUG_LOG|DEBUG_PRINT;
	clearDeb();
	const char* rFile="norm~raw.r";
	const char* R_cmd="\"C:/Program Files/R/R-3.6.0/bin/x64/Rscript.exe\"";
	const char* resP="../res";
	char b[3000], cwd[3000];
	getcwd(cwd, sizeof(cwd));

	chdir(resP);
	sprintf(b,"%s %s",R_cmd,rFile);
	deb(b);
	system(b);
	chdir(cwd);

	exit(0);
}

//===============================================================================================

int main(int argc, char **argv) {
//	tst();
	initSG(argc, argv);
//testIsFloat();
//testHTML();
	writeLog("====== Start ====== deb=%i\n",debugFg);


//===========================================
	Preparator();
	Correlator();
	fflush(stdout);
	writeLog("====== End ======\n");
	return 0;
}
