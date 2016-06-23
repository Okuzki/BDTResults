#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <math.h>
// For the progress bar. (unix only)
//#include "sys/ioctl.h"
//#include <unistd.h>
// ---------------------------------

#include "TString.h"
#include "TMath.h"
#include "TRandom.h"
#include "TChain.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

//const int FROM        = 501;
const int TO          = 5751;
//const int nCells      = 300;
const int NLAY        = 18;
const int NCELL       = 300;

bool test(const int hitNum, const int q, const std::vector<double>& chosen_cut, const std::vector<double>& edep, const std::vector<double>& tStart, const std::vector<int>& cellID, const std::vector<int>& layerID, const std::vector<double>& avgNeighbours, std::vector<int>& Tm);
//double testResult(const std::vector<double>& T, const std::vector<int>& hittype, const double FROM, const double TO );
double purity(const std::vector<double>& WBranch, const std::vector<int>& hittype);
double gini(const std::vector<double>& WBranch, const std::vector<int>& hittype);
double criterion(const double giniF, const double giniL, const double giniR);
double dVecSum(const std::vector<double>& S);
void checkneighbor(int map_k[NLAY][NCELL][NCELL]);


class ProgressBar {
	public:
		ProgressBar(int howManyBars);
		void setBar(int entries);
		void update(int iterator);
	private:
		////// For the progress bar. This progress bar works only in unix computers. //////
		//struct winsize wsize;
		//ioctl(STDOUT_FILENO,TIOCGWINSZ,&wsize);
		/* wsize.ws_row is the number of rows, size.ws_col is the number of columns. */
		int progressBarSize; //wsize.ws_col - 4;
		int nBars;
		int nEntry;
		int progressBarCount;
		int progress;
		//int pushingForwards = 1;
		int nSet;
		//double progressUpdate = 1.0/16.0;
		//int iftoosmallcheck = 0;
		std::stringstream buff;
		////// --------------------------------------------------------------------- //////

};

ProgressBar::ProgressBar(int howManyBars) { 	
	progressBarSize = 25; //wsize.ws_col - 4;
	nBars = 0;
	nEntry = 1;
	progressBarCount = 0;
	progress = 0;
	//int pushingForwards = 1;
	nSet = 0;
	nBars = howManyBars;
};

void ProgressBar::setBar(int entries) { 
	nSet++;
	nEntry = entries;
	if (nEntry*nBars < progressBarSize) {
		progressBarSize = nEntry*nBars;
		std::cout << nEntry << std::endl;
	}
};

void ProgressBar::update(int iterator) {
	//progressBarSize = 25;
	// ------------------------------- Progress bar -------------------------------
	if (nEntry == 0 || nBars == 0 || progressBarSize == 0) std::cout << iterator << ", " << nEntry << ", " << nBars << ", " << progressBarSize << std::endl;

	if ( iterator%((int)(nEntry*nBars)/(int)progressBarSize) == 0 ) {
		buff.str("");
		buff.clear();
		buff << "\33[2K |";  // \33[2K removes the current line.
		for (progress = 0; progress < progressBarCount; ++progress) buff << "=";
		//if ( c1 != TO - 1 ) {
		//	buff << ">";
		//	progress++;
		//}
		for ( ; progress < progressBarSize; ++progress) buff << "-";
		buff << "|\r";
		std::cout << buff.str();
		std::cout.flush();
		if (progressBarCount < progressBarSize) progressBarCount++;
	}
	// ----------------------------------------------------------------------------
}



int main(const int args, char **argv) {

	int nEvents = 0;
    int nTrees = 0;
	TString errorMessage = Form(" <nEvents> <nTrees>");
	if (args == 1) { 
        std::cerr << "Usage: " << argv[0] << errorMessage << std::endl;
        return -1;
    }
	if (args > 1) {
		std::istringstream ss(argv[1]);
		if (!(ss >> nEvents) && nEvents > -1) {
			std::cerr << "Usage: " << argv[0] << errorMessage << std::endl;
			return -1; 
		} 
	} 
	if (args > 2) {
		std::istringstream sss(argv[2]);
		if (!(sss >> nTrees) && nTrees > -1) {
			std::cerr << "Usage: " << argv[0] << errorMessage << std::endl;
			return -1; 
		}   
	} 
	if (args == 1) { 
        std::cerr << "Usage: " << argv[0] << errorMessage << std::endl;
        return -1;
    }

    const int NEVENTS = nEvents;
	const int NTREES = nTrees;
	const int FROM = nEvents + 1;

	time_t now = time(0);
   	char* dt = ctime(&now);

	TFile * file;
	TFile * sfile;

	TString BDTfilename = Form("file:/home/oskari/MyThing/Rootfiles/BDTOutputEv%dTree%d.root",NEVENTS,NTREES);
	TString signalFilename = Form("file:/home/oskari/MyThing/Rootfiles/signal.root");

	file = new TFile(BDTfilename, "READ"); //(("file:/Users/hanafi/Documents/Data/My_thing/BDT/BDTOutput.root", "READ");
	sfile = new TFile(signalFilename, "READ"); //(("file:/Users/hanafi/Documents/Data/My_thing/Rootfiles/signal.root", "READ");


	TTree * tree = (TTree*) file->Get("tree");
	TTree * stree = (TTree*) sfile->Get("tree");	

	double alphaSave = 0;
	std::vector<double> * treeStructure = 0;
	std::vector<double> * chosen_cut = 0;

	// int CdcCell_nHits = 0;
	std::vector<int> * CdcCell_layerID = 0;
	std::vector<int> * CdcCell_cellID = 0;
	std::vector<double> * CdcCell_edep = 0;
	// std::vector<double> * CdcCell_stepL = 0;
	// std::vector<double> * CdcCell_driftD = 0;
	// std::vector<double> * CdcCell_driftT = 0;
	std::vector<double> * CdcCell_tstart = 0;
	// std::vector<int> * CdcCell_posflag = 0;
	// std::vector<int> * CdcCell_nPair = 0;
	// std::vector<double> * CdcCell_t = 0;
	// std::vector<double> * CdcCell_px = 0;
	// std::vector<double> * CdcCell_py = 0;
	// std::vector<double> * CdcCell_pz = 0;
	// std::vector<double> * CdcCell_x = 0;
	// std::vector<double> * CdcCell_y = 0;
	// std::vector<double> * CdcCell_z = 0;
	// std::vector<double> * CdcCell_wx = 0;
	// std::vector<double> * CdcCell_wy = 0;
	// std::vector<double> * CdcCell_wz = 0;
	std::vector<int> * CdcCell_hittype = 0;
	// int M_nHits = 0;
	// std::vector<string> * M_volName = 0;
	// std::vector<int> * M_volID = 0;
	// std::vector<double> * M_edep = 0;
	// std::vector<double> * M_stepL = 0;
	// std::vector<double> * M_t = 0;
	// std::vector<double> * M_t = 0;
	// std::vector<double> * M_px = 0;
	// std::vector<double> * M_py = 0;
	// std::vector<double> * M_pz = 0;
	// std::vector<double> * M_x = 0;
	// std::vector<double> * M_y = 0;
	// std::vector<double> * M_z = 0;
	// std::vector<int> * M_hittype = 0;

	//std::vector<double> cells_edep(nCells,0);
	std::vector<double> edep;
	std::vector<int> layerID;
	std::vector<int> cellID;
	std::vector<int> hittype;
	std::vector<double> tStart;
	std::vector<double> T;
	std::vector<int> Tm;

	double alphaSummed = 0.0;


	tree->SetBranchAddress("treeStructure",&treeStructure);
	tree->SetBranchAddress("alphaSave",&alphaSave);
	tree->SetBranchAddress("chosen_cut",&chosen_cut);

	// stree->SetBranchAddress("CdcCell_nHits",&CdcCell_nHits);
	stree->SetBranchAddress("CdcCell_layerID",&CdcCell_layerID);
	stree->SetBranchAddress("CdcCell_cellID",&CdcCell_cellID);
	stree->SetBranchAddress("CdcCell_edep",&CdcCell_edep);
	// stree->SetBranchAddress("CdcCell_stepL",&CdcCell_stepL);
	// stree->SetBranchAddress("CdcCell_driftD",&CdcCell_driftD);
	// stree->SetBranchAddress("CdcCell_driftT",&CdcCell_driftT);
	stree->SetBranchAddress("CdcCell_tstart",&CdcCell_tstart);
	// stree->SetBranchAddress("CdcCell_posflag",&CdcCell_posflag);
	// stree->SetBranchAddress("CdcCell_nPair",&CdcCell_nPair);
	// stree->SetBranchAddress("CdcCell_t",&CdcCell_t);
	// stree->SetBranchAddress("CdcCell_px",&CdcCell_px);
	// stree->SetBranchAddress("CdcCell_py",&CdcCell_py);
	// stree->SetBranchAddress("CdcCell_pz",&CdcCell_pz);
	// stree->SetBranchAddress("CdcCell_x",&CdcCell_x);
	// stree->SetBranchAddress("CdcCell_y",&CdcCell_y);
	// stree->SetBranchAddress("CdcCell_z",&CdcCell_z);
	// stree->SetBranchAddress("CdcCell_wx",&CdcCell_wx);
	// stree->SetBranchAddress("CdcCell_wy",&CdcCell_wy);
	// stree->SetBranchAddress("CdcCell_wz",&CdcCell_wz);
	stree->SetBranchAddress("CdcCell_hittype",&CdcCell_hittype);
	// stree->SetBranchAddress("M_nHits",&M_nHits);
	// stree->SetBranchAddress("M_volName",&M_volName);
	// stree->SetBranchAddress("M_volID",&M_volID);
	// stree->SetBranchAddress("M_edep",&M_edep);
	// stree->SetBranchAddress("M_stepL",&M_stepL);
	// stree->SetBranchAddress("M_t",&M_t);
	// stree->SetBranchAddress("M_t",&M_t);
	// stree->SetBranchAddress("M_px",&M_px);
	// stree->SetBranchAddress("M_py",&M_py);
	// stree->SetBranchAddress("M_pz",&M_pz);
	// stree->SetBranchAddress("M_x",&M_x);
	// stree->SetBranchAddress("M_y",&M_y);
	// stree->SetBranchAddress("M_z",&M_z);
	// stree->SetBranchAddress("M_hittype",&M_hittype);

	int binTest = 1;
	int resolution = 200;
	TH1D *hResultReal = new TH1D("hResultReal","The blind results for background particles", binTest*resolution, -binTest, binTest);
	TH1D *hResultRealSignal  = new TH1D("hResultRealSignal","The blind results for signal particles", binTest*resolution, -binTest, binTest);
	TH2D *hEfficiencyBlind  = new TH2D("hEfficiencyBlind","Efficiency for the blind results", binTest*resolution, 0, binTest, binTest*resolution, 0, binTest);

	ProgressBar bar(3);// = new ProgressBar(2);

	int c1 = 0;
	int c2 = 0;
	int c3 = 0;

	unsigned nSignal  = 0;
	unsigned nBg    = 0;
	unsigned nRightSignalAnswer = 0;
	unsigned nWrongSignalAnswer = 0;
	unsigned nRightBgAnswer     = 0;
	unsigned nWrongBgAnswer     = 0;
	double sumNeighbour         = 0.0;
	int nNeighbour              = 0;
	std::vector<double> avgNeighbours;
	int map_k[NLAY][NCELL][NCELL];
	checkneighbor(map_k);

	double theBiggestT  = -1000.0;
	double theSmallestT = 1000.0;

	//unsigned changemeter = 0;


	
	bar.setBar(TO - FROM);

	for ( c1 = FROM; c1 < TO; ++c1) {
		stree->GetEntry(c1);

		bar.update(c1);
		/*
		// ------------------------------- Progress bar -------------------------------
		if ( c0%(((TO-FROM)*NUsed)/progressBarSize) == 0 ) {
			buff.str("");
			buff.clear();
			buff << "\33[2K |";  // \33[2K removes the current line.
			for (progress = 0; progress < progressBarCount; ++progress) buff << "=";
			//if ( c1 != TO - 1 ) {
			//	buff << ">";
			//	progress++;
			//}
			for ( ; progress < progressBarSize; ++progress) buff << "-";
			buff << "|\r";
			std::cout << buff.str();
			std::cout.flush();

			if (progressBarCount < progressBarSize) progressBarCount++;
		}
		// ----------------------------------------------------------------------------

		
		// ----------------- Progress bar -----------------
		iftoosmallcheck = (int)(progressUpdate*(TO - FROM));
		if ( iftoosmallcheck < 1 ) iftoosmallcheck = 1;
		if ( c1%iftoosmallcheck == 0 ) {
			buff.str("");
			buff.clear();
			buff << "\33[2K |";  // \33[2K removes the current line.
			for ( c2 = 0; c2 < progressBarCount; ++c2) buff << "=";
			//if ( c1 != TO - 1 ) {
			//	buff << ">";
			//	c2++;
			//}
			for ( ; c2 < 1.0/progressUpdate ; ++c2) buff << "-";
			buff << "|\r";
			std::cout << buff.str();
			std::cout.flush();

			progressBarCount++;
		}
		// ------------------------------------------------
		*/


		for ( c2 = 0; c2 < (int)(*CdcCell_edep).size(); ++c2) {
			edep.push_back((*CdcCell_edep).at(c2));
			hittype.push_back((*CdcCell_hittype).at(c2));
			layerID.push_back((*CdcCell_layerID).at(c2));
			cellID.push_back((*CdcCell_cellID).at(c2));
			tStart.push_back((*CdcCell_tstart).at(c2));

			T.push_back(0);
			Tm.push_back(0);

			// Counting the average neighbour energy deposition
			for (c3 = 0; c3 < (int)(*CdcCell_edep).size(); ++c3) {
				//if ((*CdcCell_layerID).at(c3) == 1 + (*CdcCell_layerID).at(c2) && (&map_k)[(*CdcCell_layerID).at(c2)][(*CdcCell_cellID).at(c2)][(*CdcCell_cellID).at(c3)] == 0) {
					if ((*CdcCell_layerID).at(c3) == (*CdcCell_layerID).at(c2) && ( (*CdcCell_cellID).at(c3) == (*CdcCell_cellID).at(c2) + 1 || (*CdcCell_cellID).at(c3) == (*CdcCell_cellID).at(c2) - 1 ) ) {
					sumNeighbour += (*CdcCell_edep).at(c3);
					nNeighbour++;
				}
			}
			if (nNeighbour > 0) {
				sumNeighbour /= (double)nNeighbour;
				avgNeighbours.push_back(sumNeighbour);
			}
			else {
				avgNeighbours.push_back(0.0);
			}
			sumNeighbour = 0.0;
			nNeighbour = 0;
		}



//		for ( int cellIDc = 0; cellIDc < nCells; ++cellIDc) {
//			for (c2 = 0; c2 < (*CdcCell_cellID).size(); ++c2) {
//				if ((*CdcCell_cellID).at(c2) == cellIDc) {
//					cells_edep.at(cellIDc) += (*CdcCell_edep).at(c2);
//				}
//			}
//		}

	}


	//progressBarCount = 0;

	/*
	// This was a test.
	//tree->Refresh();
	std::cerr << tree->GetEntries() << std::endl;
	for (int testiI = 25; testiI > 0; --testiI)
	{
		tree->GetEntry(testiI);
		//tree->Refresh();
		//tree->GetUpdate();
		std::cerr << tree->GetTreeNumber() << " " << testiI << "   ";
		std::cerr.flush();
	}
	std::cerr << std::endl;
	std::cerr << tree->GetEntries() << std::endl;
	*/

	bar.setBar(tree->GetEntries());

	for ( c1 = 0; c1 < tree->GetEntries(); ++c1 ) {
		tree->GetEntry(c1);

		bar.update(c1);
		/*
		// ----------------- Progress bar -----------------
		iftoosmallcheck = (int)(progressUpdate*tree->GetEntries());
		if ( iftoosmallcheck < 1 ) iftoosmallcheck = 1;
		if ( c1%iftoosmallcheck == 0 ) {
			buff.str("");
			buff.clear();
			buff << "\33[2K |";  // \33[2K removes the current line.
			for ( c2 = 0; c2 < progressBarCount; ++c2) buff << "=";
			//if ( c1 != TO - 1 ) {
			//	buff << ">";
			//	c2++;
			//}
			for ( ; c2 < 1.0/progressUpdate ; ++c2) buff << "-";
			buff << "|\r";
			std::cout << buff.str();
			std::cout.flush();

			progressBarCount++;
		}
		// ------------------------------------------------
		*/

		for ( c2 = 0; c2 < (int)edep.size(); ++c2 ) {
			if ( test( c2, (*treeStructure).at(0), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) {
				if ( test( c2, (*treeStructure).at(1), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) {
					if ( test( c2, (*treeStructure).at(3), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) {
						if ( test( c2, (*treeStructure).at(7), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) {
							if ( test( c2, (*treeStructure).at(15), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) { }
						}
						else {
							if ( test( c2, (*treeStructure).at(16), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) { }
						}
					}
					else {
						if ( test( c2, (*treeStructure).at(8), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) {
							if ( test( c2, (*treeStructure).at(17), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) { }
						}
						else {
							if ( test( c2, (*treeStructure).at(18), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) { }
						}
					}
				}
				else {
					if ( test( c2, (*treeStructure).at(4), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) {
						if ( test( c2, (*treeStructure).at(9), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) {
							if ( test( c2, (*treeStructure).at(19), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) { }
						}
						else {
							if ( test( c2, (*treeStructure).at(20), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) { }
						}
					}
					else {
						if ( test( c2, (*treeStructure).at(10), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) {
							if ( test( c2, (*treeStructure).at(21), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) { }
						}
						else {
							if ( test( c2, (*treeStructure).at(22), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) { }
						}
					}
				}
			}
			else {
				if ( test( c2, (*treeStructure).at(2), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) {
					if ( test( c2, (*treeStructure).at(5), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) {
						if ( test( c2, (*treeStructure).at(11), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) {
							if ( test( c2, (*treeStructure).at(23), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) { }
						}
						else {
							if ( test( c2, (*treeStructure).at(24), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) { }
						}
					}
					else {
						if ( test( c2, (*treeStructure).at(12), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) {
							if ( test( c2, (*treeStructure).at(25), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) { }
						}
						else {
							if ( test( c2, (*treeStructure).at(26), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) { }
						}
					}
				}
				else {
					if ( test( c2, (*treeStructure).at(6), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) {
						if ( test( c2, (*treeStructure).at(13), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) {
							if ( test( c2, (*treeStructure).at(27), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) { }
						}
						else {
							if ( test( c2, (*treeStructure).at(28), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) { }
						}
					}
					else {
						if ( test( c2, (*treeStructure).at(14), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) {
							if ( test( c2, (*treeStructure).at(29), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) { }
						}
						else {
							if ( test( c2, (*treeStructure).at(30), (*chosen_cut), edep, tStart, cellID, layerID, avgNeighbours, Tm) ) { }
						}
					}
				}
			}

			T.at(c2) += alphaSave*Tm.at(c2);
		}

		if (alphaSave > 0.0) alphaSummed += alphaSave;
		else alphaSummed -= alphaSave;
	}

	/*
	for (c2 = 0; c2 < Tm.size() -18; c2 += 19) {
			std::cout << std::setw(4) << Tm.at(c2)      << " " << std::setw(2) << hittype.at(c2)      << " ";
			std::cout << std::setw(4) << Tm.at(c2 + 1 ) << " " << std::setw(2) << hittype.at(c2 + 1 ) << " ";
			std::cout << std::setw(4) << Tm.at(c2 + 2 ) << " " << std::setw(2) << hittype.at(c2 + 2 ) << " ";
			std::cout << std::setw(4) << Tm.at(c2 + 3 ) << " " << std::setw(2) << hittype.at(c2 + 3 ) << " ";
			std::cout << std::setw(4) << Tm.at(c2 + 4 ) << " " << std::setw(2) << hittype.at(c2 + 4 ) << " ";
			std::cout << std::setw(4) << Tm.at(c2 + 5 ) << " " << std::setw(2) << hittype.at(c2 + 5 ) << " ";
			std::cout << std::setw(4) << Tm.at(c2 + 6 ) << " " << std::setw(2) << hittype.at(c2 + 6 ) << " ";
			std::cout << std::setw(4) << Tm.at(c2 + 7 ) << " " << std::setw(2) << hittype.at(c2 + 7 ) << " ";
			std::cout << std::setw(4) << Tm.at(c2 + 8 ) << " " << std::setw(2) << hittype.at(c2 + 8 ) << " ";
			std::cout << std::setw(4) << Tm.at(c2 + 9 ) << " " << std::setw(2) << hittype.at(c2 + 9 ) << " ";
			std::cout << std::setw(4) << Tm.at(c2 + 10) << " " << std::setw(2) << hittype.at(c2 + 10) << " ";
			std::cout << std::setw(4) << Tm.at(c2 + 11) << " " << std::setw(2) << hittype.at(c2 + 11) << " ";
			std::cout << std::setw(4) << Tm.at(c2 + 12) << " " << std::setw(2) << hittype.at(c2 + 12) << " ";
			std::cout << std::setw(4) << Tm.at(c2 + 13) << " " << std::setw(2) << hittype.at(c2 + 13) << " ";
			std::cout << std::setw(4) << Tm.at(c2 + 14) << " " << std::setw(2) << hittype.at(c2 + 14) << " ";
			std::cout << std::setw(4) << Tm.at(c2 + 15) << " " << std::setw(2) << hittype.at(c2 + 15) << " ";
			std::cout << std::setw(4) << Tm.at(c2 + 16) << " " << std::setw(2) << hittype.at(c2 + 16) << " ";
			std::cout << std::setw(4) << Tm.at(c2 + 17) << " " << std::setw(2) << hittype.at(c2 + 17) << " ";
			std::cout << std::setw(4) << Tm.at(c2 + 18) << " " << std::setw(2) << hittype.at(c2 + 18) << " ";
			std::cout << std::endl;
		}
	for (; c2 < Tm.size(); ++c2) {
		std::cout << std::setw(4) << Tm.at(c2) << " " << std::setw(2) << hittype.at(c2) << " ";
	}
	std::cout << std::endl;
	*/
	bar.setBar(T.size());

	for (c2 = 0; c2 < (int)T.size(); ++c2) {
		bar.update(c2);
		T.at(c2) /= alphaSummed;
		if (T.at(c2) > theBiggestT) theBiggestT = T.at(c2);
		if (T.at(c2) < theSmallestT) theSmallestT = T.at(c2);
		if (hittype.at(c2) < 1) hResultReal->Fill(T.at(c2));
		else if (hittype.at(c2) > 0) hResultRealSignal->Fill(T.at(c2));
	}


	// =============================================================================================================

	std::vector<double> WEqual((int)T.size(), 1.0/(double)T.size());
	std::vector<double> WL;
	std::vector<double> WR;
	double giniF = gini(WEqual, hittype);
	double resultCriterion = 0.0;
	double bestResultCrit = 0.0;
	double theBestCutForResults = -1.0;
	int nSignalTestCorrect = 0;
	int nBackgroundTestCorrect = 0;
	int nSignalResults = 0;
	int nBackground = 0;
	for (int cp1 = 0; cp1 < (int)hittype.size(); ++cp1) {
		if (hittype.at(cp1) > 0) nSignalResults++;
		else nBackground++;
	}
	WL.clear();
	WR.clear();

	bar.setBar((int)((theBiggestT - theSmallestT)/0.001));
	int testipaska = 0;

	for (double i = theSmallestT; i < theBiggestT; i += 0.001) {
		bar.update(testipaska++);
		for (int c = 0; c < (int)T.size(); ++c) {
			if ( T.at(c) < i ) {   // This is the important thing, this compares things. Thing thing
				WL.push_back(WEqual.at(c));
				WR.push_back(-1);
			}
			else {
				WL.push_back(-1);
				WR.push_back(WEqual.at(c));
			}
		}
		resultCriterion = criterion(giniF, gini(WL, hittype), gini(WR, hittype));
		if ( resultCriterion > bestResultCrit ) {
			bestResultCrit = resultCriterion;
			theBestCutForResults = i;
		}
		for (int j = 0; j < (int)WL.size(); ++j) {
			if (WL.at(j) > 0.0 && hittype.at(j) < 1) nBackgroundTestCorrect++;
			else if (WR.at(j) > 0.0 && hittype.at(j) > 0) nSignalTestCorrect++;
			//else std::cout << "W = 0 in testResult." << std::endl;
		}
		hEfficiencyBlind->Fill((double)nSignalTestCorrect/(double)nSignalResults, (double)nBackgroundTestCorrect/(double)nBackground);

		nSignalTestCorrect = 0;
		nBackgroundTestCorrect = 0;
		WL.clear();
		WR.clear();
	}

	std::cout << "The best cut for results: " << theBestCutForResults << std::endl;

	/*
	for (c2 = 0; c2 < (int)T.size(); ++c2) {
		if      ( T.at(c2) > theBestCutForResults ) sig.at(c2) = 1;    // This sig was not meant to be here?
		else sig.at(c2) = 0;
	}
	*/

	// =============================================================================================================




	//theBestCutForResults = testResult(T, hittype, theSmallestT, theBiggestT);

	for (c2 = 0; c2 < (int)T.size(); ++c2) {
		if (T.at(c2) > theBestCutForResults) {
			nSignal++;
			if (hittype.at(c2) > 0) nRightSignalAnswer++;
			else nWrongSignalAnswer++;
		}
		else if (T.at(c2) < theBestCutForResults) {
			nBg++;
			if (hittype.at(c2) < 1) nRightBgAnswer++;
			else nWrongBgAnswer++;
		}
	}

	if (nSignal + nBg == 0) std::cout << "No particles detected at all." << std::endl;
	else std::cout << "Total right answer ratio:      " << ((double)nRightBgAnswer+(double)nRightSignalAnswer)/((double)nSignal+(double)nBg) << std::endl;

	if (nBg == 0) std::cout << "No background particles detected at all." << std::endl;
	else std::cout << "Background right answer ratio: " << ((double)nRightBgAnswer)/(double)nBg << std::endl;

	if (nSignal == 0) std::cout << "No signal particles detected at all." << std::endl;
	else std::cout << "Signal right answer ratio:     " << ((double)nRightSignalAnswer)/(double)nSignal << std::endl;
	std::cout <<      "The best cut point:            " << theBestCutForResults << std::endl;


	TString fileName = Form("file:/home/oskari/MyThing/Rootfiles/BDTResultsOutputEv%dTree%d.root", NEVENTS, (int)tree->GetEntries());
	TFile *fout = new TFile(fileName,"RECREATE"); 
	fout->cd();

	// Write all histograms
	hResultReal->Write();
	hResultRealSignal->Write();
	hEfficiencyBlind->Write();

	fout->Close();

	std::cout << "start: " << dt;
	now = time(0);
   	dt = ctime(&now);
	std::cout << "end:   " << dt;

	return 0;

}










bool test(const int hitNum, const int q, const std::vector<double>& chosen_cut, const std::vector<double>& edep, const std::vector<double>& tStart, const std::vector<int>& cellID, const std::vector<int>& layerID, const std::vector<double>& avgNeighbours, std::vector<int>& Tm) {

	if      (q == 0) return (edep.at(hitNum) < chosen_cut.at(0) && edep.at(hitNum) > 0);                   // Edep cut check
    else if (q == 1) return (tStart.at(hitNum) < chosen_cut.at(1) && tStart.at(hitNum) > 0);               // tStart comparison
    else if (q == 2) return (avgNeighbours.at(hitNum) < chosen_cut.at(2) && avgNeighbours.at(hitNum) > 0); // Cell neighbourhood edep check
    else if (q == 3) return (layerID.at(hitNum) < chosen_cut.at(3) && layerID.at(hitNum) > -1);            // Layer cut check
	else if (q == 66) {
		Tm.at(hitNum) = 1;
		return false;
	}
	else if (q == -66) {
		Tm.at(hitNum) = -1;
		return false;
	}

	return false;

	// Maybe distance FROM one string hit to another?
}

/*
double testResult(const std::vector<double>& T, const std::vector<int>& hittype, const double FROM, const double TO ) {
	std::vector<double> WEqual((int)T.size(), 1.0/(double)T.size());
	double giniF = gini(WEqual, hittype);
	double resultCriterion = 0.0;
	double bestResultCrit = 0.0;
	double bestCutPoint = -1.0;
	int c2;
	std::vector<double> WL;
	std::vector<double> WR;

	
	//const double interval       = 0.001;
	//const int howManyTimes      = (int)((TO-FROM)/interval);
	//const double progressUpdate = (TO-FROM)/16.0;
	//double nowGoing             = FROM;
	//
	//// For progress bar
	//int progressBarCount = 0;
	//int iftoosmallcheck = 0;
	//std::stringstream buff;
	//// ----------------
	

	//std::cout << "Speder spudro." << howManyTimes << " ";
	//std::cout.flush();

	for (double i = FROM; i < TO; i += interval) {

		
		//// ----------------- Progress bar -----------------
		//iftoosmallcheck = (int)(progressUpdate*(1.0/0.001));
		//if ( iftoosmallcheck < 1 ) iftoosmallcheck = 1;
		//if ( i > nowGoing ) {
		//	nowGoing += progressUpdate;
		//	buff.str("");
		//	buff.clear();
		//	buff << "\33[2K |";  // \33[2K removes the current line.
		//	for ( c2 = 0; c2 < progressBarCount; ++c2) buff << "=";
		//	//if ( c1 != TO - 1 ) {
		//	//	buff << ">";
		//	//	c2++;
		//	//}
		//	for ( ; c2 < 15 ; ++c2) buff << "-";  //Fugit, that 15 is just arbitrarily there.
		//	buff << "|\r";
		//	std::cout << buff.str();
		//	std::cout.flush();
		//
		//	progressBarCount++;
		//}
		//// ------------------------------------------------
		

		for (int c = 0; c < (int)T.size(); ++c) {
			if ( T.at(c) < i ) {   // This is the important thing, this compares things. Thing thing
				WL.push_back(WEqual.at(c));
				WR.push_back(-1);
			}
			else {
				WL.push_back(-1);
				WR.push_back(WEqual.at(c));
			}
		}
		resultCriterion = criterion(giniF, gini(WL, hittype), gini(WR, hittype));
		if ( resultCriterion > bestResultCrit ) {
			bestResultCrit = resultCriterion;
			bestCutPoint = i;
		}
		WL.clear();
		WR.clear();
	}
	return bestCutPoint;
}
*/

double purity(const std::vector<double>& WBranch, const std::vector<int>& hittype) {
	double Ws = 0.0;
	double Wb = 0.0;
	double P = 0.0;
	double helpP = 0.0;

	for ( int c1P = 0; c1P < (int)WBranch.size(); ++c1P ) {
		helpP = WBranch.at(c1P);
		if ( helpP > 0.0 ) {
			if (hittype.at(c1P) > 0) Ws += helpP;
			else Wb += helpP;
		}
	}
	if (Ws == 0 && Wb == 0 ) P = 0.5;   // This is for safety.
	else P = Ws / (Ws + Wb);

	return P;
}

double gini(const std::vector<double>& WBranch, const std::vector<int>& hittype) {
	double P = purity(WBranch, hittype);

	return (dVecSum(WBranch)*P*(1-P));

}

double criterion(const double giniF, const double giniL, const double giniR) {
	return giniF - giniL - giniR;
}

double dVecSum(const std::vector<double>& S) {
	double result = 0.0;
	double helpS = 0.0;
	for ( int c1S = 0; c1S < (int)S.size(); ++c1S ) {
		helpS = S.at(c1S);
		if ( helpS > 0.0 ) result += helpS;
	}
	return result;
}



//
// This function is the courtesy of Chen Wu. I have changed some small details for my own needs.
//
//
//
void checkneighbor(int map_k[NLAY][NCELL][NCELL]) {
	TFile * TFile_wirepos = new TFile("file:/home/oskari/MyThing/Rootfiles/wirepos.140328.root", "READ"); //(("file:/Users/hanafi/Documents/Data/My_thing/Rootfiles/wirepos.140328.root", "READ");
	TTree * TTree_wirepos = (TTree*) TFile_wirepos->Get("t");
	const int entries = TTree_wirepos->GetEntries();
	std::vector<int> wmax(NLAY,0);
	//int map_k[NLAY][NCELL][NCELL];
	double map_xhv[NLAY][NCELL];
	double map_yhv[NLAY][NCELL];
	double map_xro[NLAY][NCELL];
	double map_yro[NLAY][NCELL];
	double errord [NLAY][NCELL];
	/*int*/double wp_wid;
	int wp_lid;
	int isSenseWire;
	double wp_xhv;
	double wp_yhv;
	double wp_xro;
	double wp_yro;
	TTree_wirepos->SetBranchAddress("isSenseWire",&isSenseWire);
	TTree_wirepos->SetBranchAddress("LayerID",&wp_lid);
	TTree_wirepos->SetBranchAddress("CellID",&wp_wid);
	TTree_wirepos->SetBranchAddress("xd",&wp_xro);
	TTree_wirepos->SetBranchAddress("yd",&wp_yro);
	TTree_wirepos->SetBranchAddress("xu",&wp_xhv);
	TTree_wirepos->SetBranchAddress("yu",&wp_yhv);
	for (int i = 0; i<NLAY; i++){
		wmax.at(i) = -1;
	}

	for (int i = 0; i<entries; i++){
		TTree_wirepos->GetEntry(i);
		if (!isSenseWire) continue;
		if (wp_lid>=1&&wp_lid<=17){
			map_xhv[wp_lid][(int)wp_wid] = wp_xhv;
			map_yhv[wp_lid][(int)wp_wid] = wp_yhv;
			map_xro[wp_lid][(int)wp_wid] = wp_xro;
			map_yro[wp_lid][(int)wp_wid] = wp_yro;
			errord[wp_lid][(int)wp_wid] = 0.2;
			if (wmax.at(wp_lid)<wp_wid) wmax.at(wp_lid) = wp_wid;
		}
	}
	for (int k = 1; k<NLAY-1; k++){
		//printf("layer %d & %d\n",k,k+1);
		for (int i = 0; i<=wmax.at(k); i++){
			for (int j = 0; j<=wmax.at(k+1); j++){
				if ((-(map_xro[k+1][j]-map_xhv[k+1][j])+(map_xro[k][i]-map_xhv[k][i]))){
					map_k[k][i][j] = ((map_xro[k+1][j]+map_xhv[k+1][j])-(map_xro[k][i]+map_xhv[k][i]))/(-(map_xro[k+1][j]-map_xhv[k+1][j])+(map_xro[k][i]-map_xhv[k][i]));
				}
				else{
					map_k[k][i][j] = 10;
				}
			}
		}
	}

	/*
	for (int k = 1; k<NLAY-1; k++){
		std::cout << k << ": " << std::endl;
		for (int i = 0; i<=wmax.at(k); i++){
			std::cout << i << ": ";
			for (int j = 0; j<=wmax[k+1]; j++){
				std::cout << map_k.at(k).at(i).at(j) << " ";
			}
			std::cout << std::endl;
			std::cout << std::endl;
		}
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
	}
	*/

	//return ***map_k;

}


