// This file is under public domain.
#ifndef DFTBUILDER_HPP
#define DFTBUILDER_HPP



//these are dependences, keep them
#include "commoninclude.hpp"
#include <mkl_types.h>

#include <armadillo>
#include "sparsemat.hpp"
typedef std::pair<iint,iint> At; // also defined in commoninclude.hpp
#include "mklspmat.hpp"

#include "syspara.hpp"
#include "lattice.h"
#include "electromag.hpp"
#include "excorr.hpp"





/*find Ef and occupations given total particle number, the only formula is fermi-dirac distribution, the return is found or not, depending on particle number conservation*/
bool findOccupationNumber(int nee_want, int nee_get, const double* ees, double& eehigh, double& eelow, double& kBt,double& ef, double* occu, double& deltaN,
	int& iCount) {
	using namespace std;
	bool FERMI_ENERGY_STATUS=false;
	ef = ( ees[nee_want-1] + ees[nee_want] ) / 2.0;//Firstly try one Fermi Energy, then try others in the vincinity
	deltaN=0.0;
	for (int it = 0; it < nee_get; ++it) { 
		occu[it] = 1. / (1 + exp((ees[it] - ef) / kBt)); deltaN += occu[it];
	}
	deltaN-=nee_want;
	double deltaEnergy = 5.0*abs(eehigh - eelow)+100.0*kBt;//safety comes first, (eehigh - eelow) must be a large number
	iCount = 0;
	do {
		iCount += 1;
		if (iCount > 500) {		
			if (abs(deltaN) < deltaElectronNumToleranceReal) {//Try to find the most accurate occu possible, but also accept if a flexible bar is met
				FERMI_ENERGY_STATUS = true;
				//cout << "deltaN:" << deltaN << endl;
				//cout << "iCount:" << iCount << endl;
			}
			else {
				//cout << "deltaNError:" << deltaN << endl;
				//cout << "iCountError:" << iCount << endl;//happens when fail, should do something here
			}
			return FERMI_ENERGY_STATUS;
		}
		deltaEnergy = deltaEnergy * 0.5;
		if (deltaN>0){
		ef = ef - deltaEnergy;
		}else if(deltaN < 0) {
		ef = ef + deltaEnergy;
		}

		deltaN = 0.0;
		for (int it = 0; it < nee_get; ++it) {
			occu[it] = 1. / (1 + exp((ees[it] - ef) / kBt)); deltaN += occu[it];
		}
		deltaN -= nee_want;

	} while ( abs(deltaN) > deltaElectronNumTolerance || iCount <= 5 );

	FERMI_ENERGY_STATUS = true;
	return FERMI_ENERGY_STATUS;//happens when success
}





/*yayun 240130, find Ef and occupations given total particle number, the only formula is fermi-dirac distribution, the return is found or not, depending on particle number conservation*/
bool findOccupationNumberSorted(int nee_want, int nee_get, arma::rowvec ees, double& eehigh, double& eelow, double& kBt,double& ef, double* occu, double& deltaN,
	int& iCount) {
	using namespace std;
	bool FERMI_ENERGY_STATUS=false;
	ef = ( ees[nee_want-1] + ees[nee_want] ) / 2.0;//Firstly try one Fermi Energy, then try others in the vincinity
	deltaN=0.0;
	for (int it = 0; it < nee_get; ++it) { 
		occu[it] = 1. / (1 + exp((ees[it] - ef) / kBt)); deltaN += occu[it];
	}
	deltaN-=nee_want;
	double deltaEnergy = 5.0*abs(eehigh - eelow)+100.0*kBt;//safety comes first, (eehigh - eelow) must be a large number
	iCount = 0;
	do {
		iCount += 1;
		if (iCount > 500) {		
			if (abs(deltaN) < deltaElectronNumToleranceReal) {//Try to find the most accurate occu possible, but also accept if a flexible bar is met
				FERMI_ENERGY_STATUS = true;
				//cout << "deltaN:" << deltaN << endl;
				//cout << "iCount:" << iCount << endl;
			}
			else {
				//cout << "deltaNError:" << deltaN << endl;
				//cout << "iCountError:" << iCount << endl;//happens when fail, should do something here
			}
			return FERMI_ENERGY_STATUS;
		}
		deltaEnergy = deltaEnergy * 0.5;
		if (deltaN>0){
		ef = ef - deltaEnergy;
		}else if(deltaN < 0) {
		ef = ef + deltaEnergy;
		}

		deltaN = 0.0;
		for (int it = 0; it < nee_get; ++it) {
			occu[it] = 1. / (1 + exp((ees[it] - ef) / kBt)); deltaN += occu[it];
		}
		deltaN -= nee_want;

	} while ( abs(deltaN) > deltaElectronNumTolerance || iCount <= 5 );

	FERMI_ENERGY_STATUS = true;
	return FERMI_ENERGY_STATUS;//happens when success
}




/*find Ef and occupations and make sure only one orbital has c<=Cmin*/
int findOccupationNumberCminLocate(int nee_want, int nee_get, const double* ees, double& eehigh, double& eelow, double& kBt, double& ef, double* occu, double& deltaN, int&iCount, int& nee_fix) {
	using namespace std;
	//using namespace MAIN; wrong
	int FERMI_ENERGY_STATUS_NUM = 0;
	//mainout << "test array2:" << ees[0] << "---" << occu[0] << endl;
	if (findOccupationNumber(nee_want, nee_get, ees, eehigh, eelow, kBt, ef, occu, deltaN, iCount)==false) {
		//cout << "test array3:" << ees[0] << "---" << occu[0] << endl;
		eelow = ees[0] - (ees[nee_get - 1] - ees[0]);
		eehigh = ees[nee_get - 1] + CyclotronGap;//not sure
        return FERMI_ENERGY_STATUS_NUM=1;//false means fail to find fermi etc, some accident happens: have not seen this before
	}

	//select orbitals to allow only one or none c<~Cmin


    /*upper boundary safe distance can be set by condisering orbital number, as well as nee_get */
	double val = ef + occuMinimumEkBTFirst * kBt;/*This is used to determine an energy value where occupation is eqaul to a threshold Cmin, which exists in my imagination only;
	also is the lower limit for eehigh, additionally higher energies are added to make sure there will be ees higher than val in the next iteration*/
	double valTolerance = 0.001 * kBt;//This is to find ees with energy >= val, but use ees>val-valTolerance in case that energy >= val is missed due to numerical error


	if (eehigh <= val)/* this means energy range larger than threshold, so maybe not enough orbitals */ {
		eelow = ees[0] - (ees[nee_want - 1] - ees[0]);//minus something like a "band width" that contains the number of orbitals
		eehigh = ef + occuMinimumEkBTSecond * kBt; // go to the energy upper bounds calculated by current Ef, this must be enough since rediagonalize the same Hamiltonian
		return FERMI_ENERGY_STATUS_NUM =2;// need re-diagonalize
	}
	else {
		int CminNumber = -1;
		for (int it = 0; it < nee_get; ++it) {
			if (ees[it] >= val - valTolerance) { CminNumber = it; break; }//find the first/lowest energy orbital to have occupation less than Cmin
		}
		if (CminNumber == -1) {//last orbital has c>=Cmin, but it is accidentally good to consider all orbitals in occu calculation
			eelow = ees[0] - (ees[nee_want - 1] - ees[0]);//minus something like a "band width" that contains the number of orbitals
			eehigh = ef + occuMinimumEkBTSecond * kBt; // go to the energy upper bounds calculated by current Ef, guess should be enough
			nee_fix = nee_get;
			return FERMI_ENERGY_STATUS_NUM = 3;//go to next iteration
		}
		else {
			nee_fix = CminNumber + 1;//CminNumber+1 is now the fixed total number of orbitals
			if (findOccupationNumber(nee_want, nee_fix, ees, eehigh, eelow, kBt, ef, occu, deltaN, iCount)) {
				eelow = ees[0] - (ees[nee_want - 1] - ees[0]);//minus something like a "band width" that contains the number of orbitals
				eehigh = ef + occuMinimumEkBTSecond * kBt; // go to the energy upper bounds calculated by current Ef, guess should be enough
				return FERMI_ENERGY_STATUS_NUM = 4;//go to next iteration
			}
		}
    }
	return FERMI_ENERGY_STATUS_NUM;//should not occur
}



/*find Ef and occupations and make sure only one orbital has c<=Cmin*/
int findOccuCminLocateOneaboveEf(int nee_want, int nee_get, const double* ees, double& eehigh, double& eelow, double& kBt, double& ef, double* occu, double& deltaN, int& iCount, int& nee_fix,
	const double& PotentialMin) {
	using namespace std;
	//using namespace MAIN;
	int FERMI_ENERGY_STATUS_NUM = 0;
	//mainout << "test array2:" << ees[0] << "---" << occu[0] << endl;
	if (findOccupationNumber(nee_want, nee_get, ees, eehigh, eelow, kBt, ef, occu, deltaN, iCount) == false) {
		mainout << "nee_get: " << nee_get << " nee_want: " << nee_want << "---" << endl;
		mainout << " occu[0]: " << occu[0] <<" occu[nee_get - 1]: " << occu[nee_get - 1] << endl;
		mainout << " ees[0]: " << ees[0] <<" ees[nee_get - 1]: " << ees[nee_get - 1] << endl;
		eelow = max(ees[0] - (ees[nee_get - 1] - ees[0]), PotentialMin);
		if (nee_get>1){
			eehigh = ees[nee_get - 1] + (ees[nee_get - 1] - ees[0])*(nee_want * 2.0 - nee_get)/(nee_get + nee_want);
		}else{
			eehigh = ees[nee_get - 1] + CyclotronGap*5;
		}
		/*not sure, 
		both the following problems from using continueFromfile
		yayun,231230, here max of (ees[nee_get - 1] - ees[0])) ensures a double growth, 2^10=1024 should be good enough
		yayun,231231, add a coefficient *(nee_want * 1.1 - nee_get)/nee_get*3.0 because met error: FEAST error: 3	Warning	Size of the subspace m0 is too small (m0<m)
		*/
		mainout << "false means fail to find fermi etc, it can happen if # of eigenstates is smaller than nee_want" << endl;
		return FERMI_ENERGY_STATUS_NUM = 1;//false means fail to find fermi etc, some accident happens: have not seen this before;
		//see this error message 231230, see file F:\sharedata01\yhu\workspace\data\CFDFTsquare\tri_231230. very weird
	}

	//select orbitals to allow only one or none c<~Cmin


	/*yayun, 20240126, to solve problem of "Assertion `Nsubspace >= Np' failed" which happens many times 
	see E:\sharedata01\yhu\workspace\data\CFDFTsquare\LXfracFixNee240116\30
	I here try to set up a minimum temperature  for the guess of the engenvalue range, because I am afraid that 
	it is the too small temperature of the system that causes the problems, for example suddenly miss one orbital due to small temerature occuMinimumEkBTFirst * kBt
	
	*/
	double kBt_minimum = 0.001;//yayun, 20240126, set minimum
	double kBt_guessRange = max(kBt_minimum, kBt);//yayun, 20240128, used here below, was wrong on 20240126

	/*upper boundary safe distance can be set by condisering orbital number, as well as nee_get */
	double val = ef + occuMinimumEkBTFirst * kBt;/*This is used to determine an energy value where occupation is eqaul to a threshold Cmin, which exists in my imagination only;
	also is the lower limit for eehigh, additionally higher energies are added to make sure there will be ees higher than val in the next iteration*/
	double valTolerance = 0.001 * kBt;//This is to find ees with energy >= val, but use ees>val-valTolerance in case that energy >= val is missed due to numerical error


	if (eehigh <= val)/* this means energy range larger than threshold, so maybe not enough orbitals */ {
		/*first deal with the unknown upper eigenstates, be careful with safe distance; if this (usualy) happens when a Ef in the (large) gap, then not worry much
		about the physics; A note: maybe fixNum should always keep one more extra particle*/
		if (nee_want == nee_get && eehigh > ees[nee_want - 1] + 2.0 * occuMinimumEkBTFirst * kBt/*make sure also have enough safe distance, 
            in case some upper orbitals are missed; extra higher-energy upper have less than Cmin; this condition is the main concern in the scenario of nee_want == nee_get*/
			&& /*last occupation close to 1 and Ef a lot larger; Ef will naturally go very large since require total number conservation; 
			   this is almost equivalent to the first condition of nee_want == nee_get*/
            (occu[nee_want - 1] > 1.0 - occuMinimumSecond || ees[nee_want - 1] < ef - 2.0*occuMinimumEkBTSecond * kBt 
	        || 1/*these two occu or ees conditions not necessary and even wrong, 210411 change*/)){
			//if this happens, fermi energy can be too high, should avoid by using a reasonably large fermi energy
			eelow = ees[0] - (ees[nee_want - 1] - ees[0]);//minus something like a "band width" that contains the number of orbitals
			//eehigh = ees[nee_want - 1] + 2.5 * occuMinimumEkBTSecond * kBt; 
			eehigh = ees[nee_want - 1] + 2.5 * occuMinimumEkBTSecond * kBt_guessRange; 
			/*go to the energy upper bounds calculated by last eigen plus enough safe distance, this ensures (assuming Ef stay in the middle of the last and the first upper)
			that even if there is another upper orbital after rediagonalize the same Hamiltonian, its occupation is less than Cmin*/
			ef= ees[nee_want - 1] + 2.0 * occuMinimumEkBTSecond * kBt;/*some funny large Ef occurs because nee_get=nee_want, Ef can go crazily large; artifically reset Ef*/
			for (int it = 0; it < nee_get; ++it) {
				occu[it] = 1.0; /*artifically reset occupation*/
			}
			deltaN = 0.0;
			nee_fix = nee_get;/*artifically reset everything that has to do with occupation*/
			mainout << "nee_want == nee_get happens and this is a neutral thing. " << endl;
			return FERMI_ENERGY_STATUS_NUM = 5;//go to next iteration
		}else{
		   eelow = ees[0] - (ees[nee_want - 1] - ees[0]);//minus something like a "band width" that contains the number of orbitals
		   //eehigh = ef + occuMinimumEkBTSecond * kBt; // go to the energy upper bounds calculated by current Ef, this must be enough since rediagonalize the same Hamiltonian
		   eehigh = ef + occuMinimumEkBTSecond * kBt_guessRange;
		   return FERMI_ENERGY_STATUS_NUM = 2;// need re-diagonalize
		}
	}
	else {
		int CminNumber = -1;
		for (int it = 0; it < nee_get; ++it) {
			if (ees[it] >= val - valTolerance) { CminNumber = it; break; }//find the first/lowest energy orbital to have occupation less than Cmin
		}
		if (CminNumber == -1) {//last orbital has c>=Cmin, but it is accidentally good to consider all orbitals in occu calculation
			eelow = ees[0] - (ees[nee_want - 1] - ees[0]);//minus something like a "band width" that contains the number of orbitals
			//eehigh = ef + occuMinimumEkBTSecond * kBt; // go to the energy upper bounds calculated by current Ef, guess should be enough
			eehigh = ef + occuMinimumEkBTSecond * kBt_guessRange;
			nee_fix = nee_get;
			mainout << "FERMI_ENERGY_STATUS_NUM = 3" << endl;
			return FERMI_ENERGY_STATUS_NUM = 3;//go to next iteration
		}
		else {
			nee_fix = CminNumber + 1;//CminNumber+1 is now the fixed total number of orbitals
			if (findOccupationNumber(nee_want, nee_fix, ees, eehigh, eelow, kBt, ef, occu, deltaN, iCount)) {
				eelow = ees[0] - (ees[nee_want - 1] - ees[0]);//minus something like a "band width" that contains the number of orbitals
				//eehigh = ef + occuMinimumEkBTSecond * kBt; // go to the energy upper bounds calculated by current Ef, guess should be enough
				eehigh = ef + occuMinimumEkBTSecond * kBt_guessRange;
				/*another candidate is the known and dropped eigenvalue: ees[nee_get-1], which is beyond nee_fix, 
				and presumbably beyond nee_fix again in the next iteration, so can be used as reference for upper bound of eigenvalue range
				To ensure the safety distance, add a 3-particle extra energy predicted from below 
				*/
				eehigh = min(ees[nee_get-1] + (ees[nee_get-1] - ees[nee_get-4]), eehigh);
				mainout << "FERMI_ENERGY_STATUS_NUM = 4" << endl;
				return FERMI_ENERGY_STATUS_NUM = 4;//go to next iteration
			}
		}
	}
	return FERMI_ENERGY_STATUS_NUM;//should not occur
}



//yayun, 20240127
/*find Ef and occupations and make sure only one orbital has c<=Cmin*/
int findOccuCminLocateOneaboveEf_tellExpect(int& nee_expect_temp, int nee_want, int nee_get, const double* ees, double& eehigh, double& eelow, double& kBt, double& ef, double* occu, double& deltaN, int& iCount, int& nee_fix,
	const double& PotentialMin) {
	using namespace std;
	//using namespace MAIN;
	int FERMI_ENERGY_STATUS_NUM = 0;
	//mainout << "test array2:" << ees[0] << "---" << occu[0] << endl;
	if (findOccupationNumber(nee_want, nee_get, ees, eehigh, eelow, kBt, ef, occu, deltaN, iCount) == false) {
		mainout << "nee_get: " << nee_get << " nee_want: " << nee_want << "---" << endl;
		mainout << " occu[0]: " << occu[0] <<" occu[nee_get - 1]: " << occu[nee_get - 1] << endl;
		mainout << " ees[0]: " << ees[0] <<" ees[nee_get - 1]: " << ees[nee_get - 1] << endl;
		eelow = max(ees[0] - (ees[nee_get - 1] - ees[0]), PotentialMin);
		//if this enters, nee_want < nee_get should be satisfied
		if (nee_get>1){
			if (nee_want - nee_get < 10 && (nee_want - nee_get)/nee_want < 0.1 && nee_get > 50){// if nee_get and nee_want almost the same, and system is big
				eehigh = ees[nee_get - 1] + (ees[nee_get - 1] - ees[nee_get - 20]);
			}else{
				eehigh = ees[nee_get - 1] + (ees[nee_get - 1] - ees[0])*(nee_want * 2.0 - nee_get)/(nee_get + nee_want);
			}
		}else{
			eehigh = ees[nee_get - 1] + CyclotronGap*5;
		}
		/*not sure, 
		both the following problems from using continueFromfile
		yayun,231230, here max of (ees[nee_get - 1] - ees[0])) ensures a double growth, 2^10=1024 should be good enough
		yayun,231231, add a coefficient *(nee_want * 1.1 - nee_get)/nee_get*3.0 because met error: FEAST error: 3	Warning	Size of the subspace m0 is too small (m0<m)
		*/
		mainout << "false means fail to find fermi etc, it can happen if # of eigenstates is smaller than nee_want" << endl;
		return FERMI_ENERGY_STATUS_NUM = 1;//false means fail to find fermi etc, some accident happens: have not seen this before;
		//see this error message 231230, see file F:\sharedata01\yhu\workspace\data\CFDFTsquare\tri_231230. very weird
	}

	//select orbitals to allow only one or none c<~Cmin


	/*yayun, 20240126, to solve problem of "Assertion `Nsubspace >= Np' failed" which happens many times 
	see E:\sharedata01\yhu\workspace\data\CFDFTsquare\LXfracFixNee240116\30
	I here try to set up a minimum temperature  for the guess of the engenvalue range, because I am afraid that 
	it is the too small temperature of the system that causes the problems, for example suddenly miss one orbital due to small temerature occuMinimumEkBTFirst * kBt
	
	*/
	double kBt_minimum = 0.001;//yayun, 20240126, set minimum
	double kBt_guessRange = max(kBt_minimum, kBt);//yayun, 20240128, used here below, was wrong on 20240126

	/*upper boundary safe distance can be set by condisering orbital number, as well as nee_get */
	double val = ef + occuMinimumEkBTFirst * kBt;/*This is used to determine an energy value where occupation is eqaul to a threshold Cmin, which exists in my imagination only;
	also is the lower limit for eehigh, additionally higher energies are added to make sure there will be ees higher than val in the next iteration*/
	double valTolerance = 0.001 * kBt;//This is to find ees with energy >= val, but use ees>val-valTolerance in case that energy >= val is missed due to numerical error


	if (eehigh <= val)/* this means energy range larger than threshold, so maybe not enough orbitals */ {
		/*first deal with the unknown upper eigenstates, be careful with safe distance; if this (usualy) happens when a Ef in the (large) gap, then not worry much
		about the physics; A note: maybe fixNum should always keep one more extra particle*/
		if (nee_want == nee_get && eehigh > ees[nee_want - 1] + 2.0 * occuMinimumEkBTFirst * kBt/*make sure also have enough safe distance, 
            in case some upper orbitals are missed; extra higher-energy upper have less than Cmin; this condition is the main concern in the scenario of nee_want == nee_get*/
			&& /*last occupation close to 1 and Ef a lot larger; Ef will naturally go very large since require total number conservation; 
			   this is almost equivalent to the first condition of nee_want == nee_get*/
            (occu[nee_want - 1] > 1.0 - occuMinimumSecond || ees[nee_want - 1] < ef - 2.0*occuMinimumEkBTSecond * kBt 
	        || 1/*these two occu or ees conditions not necessary and even wrong, 210411 change*/)){
			//if this happens, fermi energy can be too high, should avoid by using a reasonably large fermi energy
			eelow = ees[0] - (ees[nee_want - 1] - ees[0]);//minus something like a "band width" that contains the number of orbitals
			//eehigh = ees[nee_want - 1] + 2.5 * occuMinimumEkBTSecond * kBt; 
			eehigh = ees[nee_want - 1] + 2.5 * occuMinimumEkBTSecond * kBt_guessRange; 
			/*go to the energy upper bounds calculated by last eigen plus enough safe distance, this ensures (assuming Ef stay in the middle of the last and the first upper)
			that even if there is another upper orbital after rediagonalize the same Hamiltonian, its occupation is less than Cmin*/

			//yayun, 20240127, reverse sign
			if ( eehigh < ees[nee_get-1] + (ees[nee_get-1] - ees[nee_get-20]) ){//dos in upper temerature bound known and far
				eehigh = ees[nee_get-1] + (ees[nee_get-1] - ees[nee_get-20]);//notice that sign reverses here, prefer more orbitals	when not sure				
			}		
			nee_expect_temp = nee_get + 100;//yayun, 20240127, it happens that Nsubspace is too small when FERMI_ENERGY_STATUS_NUM == 4 "E:\sharedata01\yhu\workspace\data\CFDFTsquare\LXfracFixNee240116\28\slurm-35187_1.out"	

			ef= ees[nee_want - 1] + 2.0 * occuMinimumEkBTSecond * kBt;/*some funny large Ef occurs because nee_get=nee_want, Ef can go crazily large; artifically reset Ef*/
			for (int it = 0; it < nee_get; ++it) {
				occu[it] = 1.0; /*artifically reset occupation*/
			}
			deltaN = 0.0;
			nee_fix = nee_get;/*artifically reset everything that has to do with occupation*/
			mainout << "nee_want == nee_get happens and this is a neutral thing. " << endl;
			return FERMI_ENERGY_STATUS_NUM = 5;//go to next iteration
		}else{
		   eelow = ees[0] - (ees[nee_want - 1] - ees[0]);//minus something like a "band width" that contains the number of orbitals
		   //eehigh = ef + occuMinimumEkBTSecond * kBt; // go to the energy upper bounds calculated by current Ef, this must be enough since rediagonalize the same Hamiltonian
		   eehigh = ef + occuMinimumEkBTSecond * kBt_guessRange;
		   return FERMI_ENERGY_STATUS_NUM = 2;// need re-diagonalize
		}
	}
	else {
		int CminNumber = -1;
		for (int it = 0; it < nee_get; ++it) {
			if (ees[it] >= val - valTolerance) { CminNumber = it; break; }//find the first/lowest energy orbital to have occupation less than Cmin
		}
		if (CminNumber == -1) {//last orbital has c>=Cmin, but it is accidentally good to consider all orbitals in occu calculation
			eelow = ees[0] - (ees[nee_want - 1] - ees[0]);//minus something like a "band width" that contains the number of orbitals
			//eehigh = ef + occuMinimumEkBTSecond * kBt; // go to the energy upper bounds calculated by current Ef, guess should be enough
			eehigh = ef + occuMinimumEkBTSecond * kBt_guessRange;
			nee_fix = nee_get;
			//yayun, 20240127
			if ( eehigh > ees[nee_get-1] + (ees[nee_get-1] - ees[nee_get-20]) ){//dos in upper temerature bound known and far
				eehigh = ees[nee_get-1] + (ees[nee_get-1] - ees[nee_get-20]);					
			}		
			nee_expect_temp = nee_get + 100;//yayun, 20240127, it happens that Nsubspace is too small when FERMI_ENERGY_STATUS_NUM == 4 "E:\sharedata01\yhu\workspace\data\CFDFTsquare\LXfracFixNee240116\28\slurm-35187_1.out"	

			mainout << "FERMI_ENERGY_STATUS_NUM = 3" << endl;
			return FERMI_ENERGY_STATUS_NUM = 3;//go to next iteration
		}
		else {
			nee_fix = CminNumber + 1;//CminNumber+1 is now the fixed total number of orbitals
			if (findOccupationNumber(nee_want, nee_fix, ees, eehigh, eelow, kBt, ef, occu, deltaN, iCount)) {
				eelow = ees[0] - (ees[nee_want - 1] - ees[0]);//minus something like a "band width" that contains the number of orbitals
				//eehigh = ef + occuMinimumEkBTSecond * kBt; // go to the energy upper bounds calculated by current Ef, guess should be enough
				eehigh = ef + occuMinimumEkBTSecond * kBt_guessRange;
				/*another candidate is the known and dropped eigenvalue: ees[nee_get-1], which is beyond nee_fix, 
				and presumbably beyond nee_fix again in the next iteration, so can be used as reference for upper bound of eigenvalue range
				To ensure the safety distance, add a 3-particle extra energy predicted from below 
				*/
				if ( eehigh > ees[nee_fix-1] + (ees[nee_fix-1] - ees[nee_fix-20]) ){//dos in upper temerature bound known and far
					eehigh = ees[nee_fix-1] + (ees[nee_fix-1] - ees[nee_fix-20]);					
				}		
				nee_expect_temp = nee_get + 100;//yayun, 20240127, it happens that Nsubspace is too small when FERMI_ENERGY_STATUS_NUM == 4 "E:\sharedata01\yhu\workspace\data\CFDFTsquare\LXfracFixNee240116\28\slurm-35187_1.out"	
				mainout << "FERMI_ENERGY_STATUS_NUM = 4" << endl;
				return FERMI_ENERGY_STATUS_NUM = 4;//go to next iteration
			}
		}
	}
	return FERMI_ENERGY_STATUS_NUM;//should not occur
}



//yayun, 20240130
/*find Ef and occupations and make sure only one orbital has c<=Cmin*/
int findOccukBtGuessRange_Sorted(int& nee_expect_temp, int nee_want, int nee_get, arma::rowvec ees, double& eehigh, double& eelow, double& kBt, double& ef, double* occu, double& deltaN, int& iCount, int& nee_fix,
	const double PotentialMin) {
	using namespace std;
	//using namespace MAIN;
	int FERMI_ENERGY_STATUS_NUM = 0;
	//mainout << "test array2:" << ees[0] << "---" << occu[0] << endl;
	if (findOccupationNumberSorted(nee_want, nee_get, ees, eehigh, eelow, kBt, ef, occu, deltaN, iCount) == false) {
		mainout << "nee_get: " << nee_get << " nee_want: " << nee_want << "---" << endl;
		mainout << " occu[0]: " << occu[0] <<" occu[nee_get - 1]: " << occu[nee_get - 1] << endl;
		mainout << " ees[0]: " << ees[0] <<" ees[nee_get - 1]: " << ees[nee_get - 1] << endl;
		eelow = max(eelow - (ees[nee_get - 1] - ees[0]), PotentialMin);
		//if this enters, nee_want < nee_get should be satisfied
		if (nee_get>1){
			if(nee_get > 100){
				if (nee_want - nee_get < 20 && (nee_want - nee_get)/nee_want < 0.2){// if nee_get and nee_want almost the same, and system is big
					eehigh = eehigh + (ees[nee_get - 1] - ees[nee_get - 100]);//reference is eehigh itself
				}else{
					eehigh = eehigh + (eehigh - ees[0]) * (nee_want - nee_get)/nee_get * 2.0;/*not close, estimate useing dos; because I always look at energies in the gap 
					that has a relatively dos, and I should be so here since I already have >100 orbitals, thus only *2.0 should be fine. Also because not risk too many orbitals*/
				}
			}else{
				eehigh = eehigh + (eehigh - ees[0]) * max(1.0, (nee_want - nee_get)/nee_get * 2.0);//not close and few orbitals, double without much cost or risk
			}			
		}else{// no eigen found , double
			eehigh = eehigh + CyclotronGap * 5.0 + (eehigh - eelow) * 1.0;//reference is eehigh itself
		}
		/*not sure, 
		both the following problems from using continueFromfile
		yayun,231230, here max of (ees[nee_get - 1] - ees[0])) ensures a double growth, 2^10=1024 should be good enough
		yayun,231231, add a coefficient *(nee_want * 1.1 - nee_get)/nee_get*3.0 because met error: FEAST error: 3	Warning	Size of the subspace m0 is too small (m0<m)
		*/
		mainout << "false means fail to find fermi etc, it can happen if # of eigenstates is smaller than nee_want" << endl;
		return FERMI_ENERGY_STATUS_NUM = 1;//false means fail to find fermi etc, some accident happens: have not seen this before;
		//see this error message 231230, see file F:\sharedata01\yhu\workspace\data\CFDFTsquare\tri_231230. very weird
	}

	//select orbitals to allow only one or none c<~Cmin


	/*yayun, 20240126, to solve problem of "Assertion `Nsubspace >= Np' failed" which happens many times 
	see E:\sharedata01\yhu\workspace\data\CFDFTsquare\LXfracFixNee240116\30
	I here try to set up a minimum temperature  for the guess of the engenvalue range, because I am afraid that 
	it is the too small temperature of the system that causes the problems, for example suddenly miss one orbital due to small temerature occuMinimumEkBTFirst * kBt
	
	*/
	double kBt_minimum = 0.001;//yayun, 20240126, set minimum
	double kBt_guessRange = max(kBt_minimum, kBt);//yayun, 20240128, used here below, was wrong on 20240126

	/*upper boundary safe distance can be set by condisering orbital number, as well as nee_get */
	double val = ef + occuMinimumEkBTFirst * kBt;/*This is used to determine an energy value where occupation is eqaul to a threshold Cmin, which exists in my imagination only;
	also is the lower limit for eehigh, additionally higher energies are added to make sure there will be ees higher than val in the next iteration*/
	double valTolerance = 0.001 * kBt;//This is to find ees with energy >= val, but use ees>val-valTolerance in case that energy >= val is missed due to numerical error


	if (eehigh <= val)/* this means energy range larger than threshold, so maybe not enough orbitals */ {
		/*first deal with the unknown upper eigenstates, be careful with safe distance; if this (usualy) happens when a Ef in the (large) gap, then not worry much
		about the physics; A note: maybe fixNum should always keep one more extra particle*/
		if (nee_want == nee_get && eehigh > ees[nee_want - 1] + 2.0 * occuMinimumEkBTFirst * kBt/*make sure also have enough safe distance, 
            in case some upper orbitals are missed; extra higher-energy upper have less than Cmin; this condition is the main concern in the scenario of nee_want == nee_get*/
			&& /*last occupation close to 1 and Ef a lot larger; Ef will naturally go very large since require total number conservation; 
			   this is almost equivalent to the first condition of nee_want == nee_get*/
            (occu[nee_want - 1] > 1.0 - occuMinimumSecond || ees[nee_want - 1] < ef - 2.0*occuMinimumEkBTSecond * kBt 
	        || 1/*these two occu or ees conditions not necessary and even wrong, 210411 change*/)){
			//if this happens, fermi energy can be too high, should avoid by using a reasonably large fermi energy
			eelow = ees[0] - (ees[nee_want - 1] - ees[0]);//minus something like a "band width" that contains the number of orbitals
			//eehigh = ees[nee_want - 1] + 2.5 * occuMinimumEkBTSecond * kBt; 
			eehigh = ees[nee_want - 1] + 2.5 * occuMinimumEkBTSecond * kBt_guessRange; 
			/*go to the energy upper bounds calculated by last eigen plus enough safe distance, this ensures (assuming Ef stay in the middle of the last and the first upper)
			that even if there is another upper orbital after rediagonalize the same Hamiltonian, its occupation is less than Cmin*/

			//yayun, 20240127, reverse sign
			if ( eehigh < ees[nee_get-1] + (ees[nee_get-1] - ees[nee_get-20]) ){//dos in upper temerature bound known and far
				eehigh = ees[nee_get-1] + (ees[nee_get-1] - ees[nee_get-20]);//notice that sign reverses here, prefer more orbitals	when not sure				
			}		
			nee_expect_temp = nee_get + 100;//yayun, 20240127, it happens that Nsubspace is too small when FERMI_ENERGY_STATUS_NUM == 4 "E:\sharedata01\yhu\workspace\data\CFDFTsquare\LXfracFixNee240116\28\slurm-35187_1.out"	

			ef= ees[nee_want - 1] + 2.0 * occuMinimumEkBTSecond * kBt;/*some funny large Ef occurs because nee_get=nee_want, Ef can go crazily large; artifically reset Ef*/
			for (int it = 0; it < nee_get; ++it) {
				occu[it] = 1.0; /*artifically reset occupation*/
			}
			deltaN = 0.0;
			nee_fix = nee_get;/*artifically reset everything that has to do with occupation*/
			mainout << "nee_want == nee_get happens and this is a neutral thing. " << endl;
			return FERMI_ENERGY_STATUS_NUM = 5;//go to next iteration
		}else{
		   eelow = ees[0] - (ees[nee_want - 1] - ees[0]);//minus something like a "band width" that contains the number of orbitals
		   //eehigh = ef + occuMinimumEkBTSecond * kBt; // go to the energy upper bounds calculated by current Ef, this must be enough since rediagonalize the same Hamiltonian
		   eehigh = ef + occuMinimumEkBTSecond * kBt_guessRange;
		   return FERMI_ENERGY_STATUS_NUM = 2;// need re-diagonalize
		}
	}
	else {
		int CminNumber = -1;
		for (int it = 0; it < nee_get; ++it) {
			if (ees[it] >= val - valTolerance) { CminNumber = it; break; }//find the first/lowest energy orbital to have occupation less than Cmin
		}
		if (CminNumber == -1) {//last orbital has c>=Cmin, but it is accidentally good to consider all orbitals in occu calculation
			eelow = ees[0] - (ees[nee_want - 1] - ees[0]);//minus something like a "band width" that contains the number of orbitals
			//eehigh = ef + occuMinimumEkBTSecond * kBt; // go to the energy upper bounds calculated by current Ef, guess should be enough
			eehigh = ef + occuMinimumEkBTSecond * kBt_guessRange;
			nee_fix = nee_get;
			//yayun, 20240127
			if ( eehigh > ees[nee_get-1] + (ees[nee_get-1] - ees[nee_get-20]) ){//dos in upper temerature bound known and far
				eehigh = ees[nee_get-1] + (ees[nee_get-1] - ees[nee_get-20]);					
			}		
			nee_expect_temp = nee_get + 100;//yayun, 20240127, it happens that Nsubspace is too small when FERMI_ENERGY_STATUS_NUM == 4 "E:\sharedata01\yhu\workspace\data\CFDFTsquare\LXfracFixNee240116\28\slurm-35187_1.out"	

			mainout << "FERMI_ENERGY_STATUS_NUM = 3" << endl;
			return FERMI_ENERGY_STATUS_NUM = 3;//go to next iteration
		}
		else {
			nee_fix = CminNumber + 1;//CminNumber+1 is now the fixed total number of orbitals
			if (findOccupationNumberSorted(nee_want, nee_fix, ees, eehigh, eelow, kBt, ef, occu, deltaN, iCount)) {
				eelow = ees[0] - (ees[nee_want - 1] - ees[0]);//minus something like a "band width" that contains the number of orbitals
				//eehigh = ef + occuMinimumEkBTSecond * kBt; // go to the energy upper bounds calculated by current Ef, guess should be enough
				eehigh = ef + occuMinimumEkBTSecond * kBt_guessRange;
				/*another candidate is the known and dropped eigenvalue: ees[nee_get-1], which is beyond nee_fix, 
				and presumbably beyond nee_fix again in the next iteration, so can be used as reference for upper bound of eigenvalue range
				To ensure the safety distance, add a 3-particle extra energy predicted from below 
				*/
				if ( eehigh > ees[nee_fix-1] + (ees[nee_fix-1] - ees[nee_fix-20]) ){//dos in upper temerature bound known and far
					eehigh = ees[nee_fix-1] + (ees[nee_fix-1] - ees[nee_fix-20]);					
				}		
				nee_expect_temp = nee_get + 100;//yayun, 20240127, it happens that Nsubspace is too small when FERMI_ENERGY_STATUS_NUM == 4 "E:\sharedata01\yhu\workspace\data\CFDFTsquare\LXfracFixNee240116\28\slurm-35187_1.out"	
				mainout << "FERMI_ENERGY_STATUS_NUM = 4" << endl;
				return FERMI_ENERGY_STATUS_NUM = 4;//go to next iteration
			}
		}
	}
	return FERMI_ENERGY_STATUS_NUM;//should not occur
}



/*find occupations for a given Ef*/
void findOccupationNumberGivenEf(double& nee_wantdecimal, int& nee_get, const double* ees, double& kBt, double& ef, double* occu) {
	using namespace std;
	nee_wantdecimal = 0.0;
	for (int it = 0; it < nee_get; ++it) {
		occu[it] = 1. / (1 + exp((ees[it] - ef) / kBt)); 
		nee_wantdecimal += occu[it];
	}	
}

/*find occupations for a given Ef*/
void findOccupationNumberGivenEfSorted(double& nee_wantdecimal, int& nee_get, arma::rowvec ees, double& kBt, double& ef, double* occu) {
	using namespace std;
	nee_wantdecimal = 0.0;
	for (int it = 0; it < nee_get; ++it) {
		occu[it] = 1. / (1 + exp((ees[it] - ef) / kBt)); 
		nee_wantdecimal += occu[it];
	}	
}

/*find Ef and occupations and make sure only one orbital has c<=Cmin; here the Ef is fixed*/
int findOccuCminfixedEf(int& nee_want, int& nee_get, const double* ees, double& eehigh, double& eelow, double& kBt, double& ef, double* occu,
	int& nee_fix, const double& PotentialMintemp, double& lastEmin, double& nee_wantdecimal) {
	using namespace std;
	int FERMI_ENERGY_STATUS_NUM = 0;

	findOccupationNumberGivenEf(nee_wantdecimal,nee_get, ees, kBt, ef, occu);
	eelow = min( min(ees[0], lastEmin) - (ees[nee_get - 1] - ees[0]) * 0.2, PotentialMintemp);/*Probably ees[0]-bandwidth will be chosen; I assume that PotentialMintemp is the lower bound of eelow,
	such that next time I can still catch the lowest energy orbital even if I missed it by accident*/
	/*(ees[nee_get - 1] - ees[0]) * 0.2, I did not change coefficient 0.2 to 1.0 when I meet the problem of search over 10 times abortion, because the lower bound does not change in the search. Yayun 231230*/
	eehigh = ef + occuMinimumEkBTFirst * kBt;//this is convenient.
	/*when I meet the problem of search over 10 times abortion in zeroth iteration, because the upper bound does not grow quickly enough in the search. 
	I thought that the starting temperature should not be too samll, which may slow down the above update. But not really, because zeroth iteration does not come here when efmode=1. It is the efmode=1 part responsible.  Yayun 231230 */

	lastEmin = ees[0];
	nee_fix = nee_get;
	nee_want = round(nee_wantdecimal);//??
	return FERMI_ENERGY_STATUS_NUM=2;//happens when success
}


/*find Ef and occupations and make sure only one orbital has c<=Cmin; here the Ef is fixed*/
int findOccuCminfixedEfSorted(int& nee_want, int& nee_get, arma::rowvec ees, double& eehigh, double& eelow, double& kBt, double& ef, double* occu,
	int& nee_fix, const double& PotentialMintemp, double& lastEmin, double& nee_wantdecimal) {
	using namespace std;
	int FERMI_ENERGY_STATUS_NUM = 0;

	findOccupationNumberGivenEfSorted(nee_wantdecimal,nee_get, ees, kBt, ef, occu);
	eelow = min( min(ees[0], lastEmin) - (ees[nee_get - 1] - ees[0]) * 0.2, PotentialMintemp);/*Probably ees[0]-bandwidth will be chosen; I assume that PotentialMintemp is the lower bound of eelow,
	such that next time I can still catch the lowest energy orbital even if I missed it by accident*/
	/*(ees[nee_get - 1] - ees[0]) * 0.2, I did not change coefficient 0.2 to 1.0 when I meet the problem of search over 10 times abortion, because the lower bound does not change in the search. Yayun 231230*/
	eehigh = ef + occuMinimumEkBTFirst * kBt;//this is convenient.
	/*when I meet the problem of search over 10 times abortion in zeroth iteration, because the upper bound does not grow quickly enough in the search. 
	I thought that the starting temperature should not be too samll, which may slow down the above update. But not really, because zeroth iteration does not come here when efmode=1. It is the efmode=1 part responsible.  Yayun 231230 */

	lastEmin = ees[0];
	nee_fix = nee_get;
	nee_want = round(nee_wantdecimal);//??
	return FERMI_ENERGY_STATUS_NUM=2;//happens when success
}



bool shiftEboudFirst(int nee_want, int& nee_get, const double* ees, double& eehigh, double& eelow, double& kBt, double& ef, double* occu, const bool& abool, bool& bbool){
	bool SHIFT_EBOUND_STATUS = false;
	//occuMinimumFirst;
	//occuMinimumEkBTFirst;
	//occuMinimumSecond;
	//occuMinimumEkBTSecond;
	//CyclotronGap;

	/*
	three key concepts:
	a "band width" defined by the fixed number orbitals contained, e.g., obtained by setbackwardNumber
	an energy distance equialence of decay of occupation number, e.g., occuMinimumEkBTFirst * kBt 
	also the CyclotronGap
	*/
	
	/*low boundary safe distance can be set by condisering band width from orbital number, temeprature, and energy gap */
	int setbackwardNumber = (nee_want< nee_get)? nee_want : nee_get;//pick the small one, usually nee_want is smaller; this is something like a "band width" that contains certain number of orbitals 
	double eLboundNumber = ees[setbackwardNumber - 1] - ees[0];
	double eLboundTemperature = occuMinimumEkBTFirst * kBt;
	double eLbound = (eLboundTemperature>eLboundNumber) ? eLboundTemperature : eLboundNumber;
	eLbound = (CyclotronGap>eLbound) ? CyclotronGap : eLbound;//pick the largest energy distance, usually eLboundNumber is larger?;
	eelow = ees[0] - eLbound;

	/*upper boundary safe distance can be set by condisering orbital number, temeprature, and energy gap, as well as nee_get */
	double val = ef + occuMinimumEkBTFirst * kBt;/*This is used to determine an energy value where occupation is eqaul to a threshold Cmin, which exists in my imagination only;
	also is the lower limit for eehigh, additionally higher energies are added to make sure there will be ees higher than val in the next iteration*/
	double valTolerance = 0.001 * kBt;//This is to find ees with energy >= val, but use ees>val-valTolerance in case that energy >= val is missed due to numerical error
	double eHbound;
	int setforwardNumber;//This is also used to calculate something like a "band width" defined by the fixed number orbitals contained
	/*Notice we can allow recalculating occupationNumber and fermiEnergy here by using only a desired number of orbitals, this is useful because occu may very sensitively depend on orbital number and 
	I do not want to calculate that many of them, so I fix it.*/
	int maxNumOrbital = 1000;
	if (nee_get > maxNumOrbital) {// here is a situation where the orbital number is too many, so smaller safe distance
		if (ees[maxNumOrbital-1] <= val)/* this means occupation of ees[999=maxNumOrbital-1] larger than threshold Cmin */ {
			eHbound = (occuMinimumEkBTSecond - occuMinimumEkBTFirst) * kBt;
			setforwardNumber = (maxNumOrbital - 100 > 0) ? maxNumOrbital - 100 : 0;
			eHbound = ((ees[maxNumOrbital - 1] - ees[setforwardNumber])>eHbound) ? (ees[maxNumOrbital - 1] - ees[setforwardNumber]) : eHbound;
			eHbound = (0.5*CyclotronGap>eHbound) ? 0.5*CyclotronGap : eHbound;
			if (abool) {				
				double deltaNdiff; int icount = -1;//this definition is only to make up for the reuiqred number of arguments for findOccupationNumber
				if (findOccupationNumber(nee_want, maxNumOrbital, ees, eehigh, eelow, kBt, ef, occu, deltaNdiff, icount)) {//This means the orbitals more than maxNumOrbital are ignored
					eehigh = ees[maxNumOrbital - 1] + eHbound;
					nee_get = maxNumOrbital;// use max to make sure smooth evolution of the relation nee_get > maxNumOrbital, such that do not miss/lose orbitals suddenly next time 
					bbool = true;
					SHIFT_EBOUND_STATUS = true;
					return SHIFT_EBOUND_STATUS;
				}
				else {
					eehigh = ees[nee_get- 1] + eHbound;// This choice happens only if above does not find fermi energy, does not seem possible; Expand search region in this case;
					return SHIFT_EBOUND_STATUS; // need re-diagonalize
				}
			}
		}
		else /* this means occupation of ees[999=maxNumOrbital-1] smaller than threshold Cmin */ {
			int miniFirstNumber = 0;
			for (int it = 0; it < nee_get; ++it) {
				if (ees[it] >= val - valTolerance) { miniFirstNumber = it; break; }//find the first/lowest energy orbital to have occupation less than Cmin
			}			
		    eHbound = (occuMinimumEkBTSecond - occuMinimumEkBTFirst) * kBt;         
			setforwardNumber = (miniFirstNumber - 50 > 0)? miniFirstNumber - 50 : 0;
			eHbound = ((ees[miniFirstNumber] - ees[setforwardNumber]) > eHbound) ? (ees[miniFirstNumber] - ees[setforwardNumber]) : eHbound;
			eHbound = (0.5*CyclotronGap>eHbound) ? 0.5*CyclotronGap : eHbound;
			eehigh = ees[miniFirstNumber] + eHbound;// use max to make sure the desired miniFirstNumber can still be found next time
			SHIFT_EBOUND_STATUS = true;
			return SHIFT_EBOUND_STATUS;
		}
	}
	else if (nee_get > 3 * nee_want && nee_get <= maxNumOrbital) {// here is a situation where the orbital number is not too many but not too few, , so moderate safe distance
		if (ees[nee_get - 2] <= val)/* this means occupation of last but one ees larger than threshold, so not enough orbitals */ {
			eHbound = ef + (occuMinimumEkBTSecond * kBt > 0.5*CyclotronGap) ? occuMinimumEkBTSecond * kBt : 0.5*CyclotronGap; // go to the energy upper bounds calculated by current Ef
			return SHIFT_EBOUND_STATUS;// need re-diagonalize
		}
		else {
			int miniFirstNumber = 0;
			for (int it = 0; it < nee_get; ++it) {
				if (ees[it] >= val - valTolerance) { miniFirstNumber = it; break; }//find the first/lowest energy orbital to have occupation less than Cmin
			}
			eHbound = (occuMinimumEkBTSecond - occuMinimumEkBTFirst) * kBt;
			setforwardNumber = (miniFirstNumber - 50 > 0) ? miniFirstNumber - 50 : 0;
			eHbound = ((ees[miniFirstNumber] - ees[setforwardNumber]) > eHbound) ? (ees[miniFirstNumber] - ees[setforwardNumber]) : eHbound;
			eHbound = (0.5*CyclotronGap > eHbound) ? 0.5*CyclotronGap : eHbound;
			eehigh = ees[miniFirstNumber] + eHbound;// use max to make sure the desired miniFirstNumber can still be found next time
			SHIFT_EBOUND_STATUS = true;
			return SHIFT_EBOUND_STATUS;
		}
	}
	else if(nee_get <= 3 * nee_want && nee_get>=1) {// here is a situation where the orbital number is too few, so larger safe distance
		if (ees[nee_get - 2] <= val)/* this means occupation of last but one ees larger than threshold, so not enough orbitals */ {
			eHbound = ef + (occuMinimumEkBTSecond * kBt > 2*CyclotronGap) ? occuMinimumEkBTSecond * kBt : 2*CyclotronGap; // go to the energy upper bounds calculated by current Ef
			return SHIFT_EBOUND_STATUS;// need re-diagonalize
		}
		else {
			int miniFirstNumber = 0;
			for (int it = 0; it < nee_get; ++it) {
				if (ees[it] >= val - valTolerance) { miniFirstNumber = it; break; }
			}
			eHbound = (occuMinimumEkBTSecond - occuMinimumEkBTFirst) * kBt;
			setforwardNumber = (miniFirstNumber - 150 > 0) ? miniFirstNumber - 150 : 0;
			eHbound = ((ees[miniFirstNumber] - ees[setforwardNumber]) > eHbound) ? (ees[miniFirstNumber] - ees[setforwardNumber]) : eHbound;
			eHbound = (2.0*CyclotronGap > eHbound) ? 2.0*CyclotronGap : eHbound;
			eehigh = ees[miniFirstNumber] + eHbound;// use max to make sure the desired miniFirstNumber can still be found next time
			SHIFT_EBOUND_STATUS = true;
			return SHIFT_EBOUND_STATUS;
		}
	}
	else if (nee_get <=0) {// This is not possible to happen given that we calculate fermiEnergy etc previously
		eHbound = ef + (occuMinimumEkBTSecond * kBt > 2 * CyclotronGap) ? occuMinimumEkBTSecond * kBt : 2 * CyclotronGap; // go to the energy upper bounds calculated by current Ef
		return SHIFT_EBOUND_STATUS;// need re-diagonalize
	}

	return SHIFT_EBOUND_STATUS;
}



/*this (1) examines if the current fermi and occu calculation is reasonable, (2) also dicedes if should re-digonalize.
(3) and make next engenvalue range guess. */
bool shiftEboudFirstOccuTemperature(int nee_want, int& nee_get, const double* ees, double& eehigh, double& eelow, double& kBt, double& ef, double* occu, const bool& abool, bool& bbool) {
	bool SHIFT_EBOUND_STATUS = false;
	/*
	three key concepts:
1.	a "band width" defined by the fixed number orbitals contained, e.g., obtained by setbackwardNumber
2.	an energy distance equialence of decay of occupation number, e.g., occuMinimumEkBTFirst * kBt
3.	also the CyclotronGap
	
	for upper bound, apply 2.kBT mostly/even only, prefer not to use 3.CyclotronGap; 
	*/
	
	/*low boundary safe distance can be set by condisering band width from orbital number, temeprature, and energy gap */
	//as as in shiftEboudFirst 
	int setbackwardNumber = (nee_want < nee_get) ? nee_want : nee_get;//pick the small one, usually nee_want is smaller; this is something like a "band width" that contains certain number of orbitals 
	double eLboundNumber = ees[setbackwardNumber - 1] - ees[0];
	double eLboundTemperature = occuMinimumEkBTFirst * kBt;
	double eLbound = (eLboundTemperature > eLboundNumber) ? eLboundTemperature : eLboundNumber;
	eLbound = (CyclotronGap > eLbound) ? CyclotronGap : eLbound;//pick the largest energy distance, usually eLboundNumber is larger?;
	eelow = ees[0] - eLbound;

	/*upper boundary safe distance can be set by condisering orbital number, temeprature, and energy gap, as well as nee_get */
	// prefer to use kBT than in shiftEboudFirst, no CyclotronGap used
	double val = ef + occuMinimumEkBTFirst * kBt;/*This is used to determine an energy value where occupation is eqaul to a threshold Cmin, which exists in my imagination only;
	also is the lower limit for eehigh, additionally higher energies are added to make sure there will be ees higher than val in the next iteration*/
	double valTolerance = 0.001 * kBt;//This is to find ees with energy >= val, but use ees>val-valTolerance in case that energy >= val is missed due to numerical error
	double eHbound;
	int setforwardNumber;//This is also used to calculate something like a "band width" defined by the fixed number orbitals contained
	/*Notice we can allow recalculating occupationNumber and fermiEnergy here by using only a desired number of orbitals, this is useful because occu may very sensitively depend on orbital number and
	I do not want to calculate that many of them, so I fix it.*/
	int maxNumOrbital = 1000;
	if (nee_get > maxNumOrbital) {// here is a situation where the orbital number is too many, so smaller safe distance
		if (ees[maxNumOrbital - 1] <= val)/* this means occupation of ees[999=maxNumOrbital-1] larger than threshold Cmin */ {
			//eHbound = (occuMinimumEkBTSecond - occuMinimumEkBTFirst) * kBt; modified from shiftEboudFirst
			setforwardNumber = (maxNumOrbital - 100 > 0) ? maxNumOrbital - 100 : 0;
			eHbound = ees[maxNumOrbital - 1] - ees[setforwardNumber];//make sure have ~100 more than maxNumOrbital
			//eHbound = ((ees[maxNumOrbital - 1] - ees[setforwardNumber]) > eHbound) ? (ees[maxNumOrbital - 1] - ees[setforwardNumber]) : eHbound;modified from shiftEboudFirst
			//eHbound = (0.5 * CyclotronGap > eHbound) ? 0.5 * CyclotronGap : eHbound; modified from shiftEboudFirst
			if (abool) {
				double deltaNdiff; int icount = -1;//this definition is only to make up for the reuiqred number of arguments for findOccupationNumber
				if (findOccupationNumber(nee_want, maxNumOrbital, ees, eehigh, eelow, kBt, ef, occu, deltaNdiff, icount)) {//This means the orbitals more than maxNumOrbital are ignored
					eehigh = ees[maxNumOrbital - 1] + eHbound;
					nee_get = maxNumOrbital;// use max to make sure smooth evolution of the relation nee_get > maxNumOrbital, such that do not miss/lose orbitals suddenly next time 
					bbool = true;
					SHIFT_EBOUND_STATUS = true;
					return SHIFT_EBOUND_STATUS;
					/*PHILOSOPHY: new ef_new will be larger than old ef_old due to less orbitals, so if val_old=ef_old+xx is larger than es[maxNumOrbital - 1], 
					then val_new=ef_new+xx >= es[maxNumOrbital - 1];
					*in next solutions, orbital numbers be ~maxNumOrbital+100, ef_next_iter_old be around last ef_new or less since ~maxNumOrbital+100> maxNumOrbital, or
					* val_old < ef_next_iter_old < val_new, this would risk ef_next_iter_old+xx >= es[maxNumOrbital - 1] and 
					as a result go to the following if-loop branch: ef_next_iter_old+xx >= es[maxNumOrbital - 1];
					which will then involve kBT to set upper limit to include even larger # orbitals beyond ~maxNumOrbital+100; 
					To walk around this problem: use the smaller upper limit that are determined from kBT and setforwardNumber in the other if-loop branch
					REASON: in the next if-loop branch, ef_next_iter_old+xx is larger than it "should be when orbital # infinite", thus eHbound_kBT = ef + occuMinimumEkBTSecond * kBt is overestimated,
					if eHbound_kBT < eHbound_setforwardNumber, then it means ees[~maxNumOrbital+100] > Cmin when the orbital # constraint of ~maxNumOrbital+100 is released.
					*/
				}
				else {
					eHbound = (1.0 * CyclotronGap > eHbound) ? 1.0 * CyclotronGap : eHbound;//modified from shiftEboudFirst
					eehigh = ees[nee_get - 1] + eHbound;// This choice happens only if above does not find fermi energy, does not seem possible; Expand search region in this case;
					return SHIFT_EBOUND_STATUS; // need re-diagonalize
				}
			}
		}
		else /* this means occupation of ees[999=maxNumOrbital-1] smaller than threshold Cmin */ {//"completely" modified from shiftEboudFirst
			eHbound = ef + occuMinimumEkBTSecond * kBt - ees[maxNumOrbital - 1];//this might be negative
			setforwardNumber = (maxNumOrbital-1 - 100 > 0) ? maxNumOrbital-1 - 100 : 0;
			eHbound = ((ees[maxNumOrbital - 1] - ees[setforwardNumber]) < eHbound) ? (ees[maxNumOrbital - 1] - ees[setforwardNumber]) : eHbound;//pick the smaller one
			//eHbound = (0.5 * CyclotronGap > eHbound) ? 0.5 * CyclotronGap : eHbound; modified from shiftEboudFirst
			eehigh = ees[maxNumOrbital - 1] + eHbound;// use min to make sure the orbital # can decrease when it "should", and can not increase when it "should not"
			SHIFT_EBOUND_STATUS = true;
			return SHIFT_EBOUND_STATUS;
		}
	}
	else if (nee_get > 3 * nee_want && nee_get <= maxNumOrbital) {// here is a situation where the orbital number is not too many but not too few, , so moderate safe distance
		if (ees[nee_get - 2] <= val)/* this means occupation of last but one ees larger than threshold Cmin, so not enough orbitals */ {
			eHbound = ef + (occuMinimumEkBTSecond * kBt > 1.0 * CyclotronGap) ? occuMinimumEkBTSecond * kBt : 1.0 * CyclotronGap; // go to the energy upper bounds calculated by current Ef
			return SHIFT_EBOUND_STATUS;// need re-diagonalize
		}
		else {
			int miniFirstNumber = 0;
			for (int it = 0; it < nee_get; ++it) {
				if (ees[it] >= val - valTolerance) { miniFirstNumber = it; break; }//find the first/lowest energy orbital to have occupation less than Cmin
			}
			eHbound = (occuMinimumEkBTSecond - occuMinimumEkBTFirst) * kBt;
			setforwardNumber = (miniFirstNumber - nee_want > 0) ? miniFirstNumber - nee_want : 0;//nee_want may be the Landau level degeneracy, such that can bear with level crossing around EF
			eHbound = ((ees[miniFirstNumber] - ees[setforwardNumber]) > eHbound) ? (ees[miniFirstNumber] - ees[setforwardNumber]) : eHbound;
			//eHbound = (0.5 * CyclotronGap > eHbound) ? 0.5 * CyclotronGap : eHbound; modified from shiftEboudFirst
			eehigh = ees[miniFirstNumber] + eHbound;// use max to make sure the desired miniFirstNumber can still be found next time
			SHIFT_EBOUND_STATUS = true;
			return SHIFT_EBOUND_STATUS;
		}
	}
	else if (nee_get <= 3 * nee_want && nee_get >= 1) {// here is a situation where the orbital number is too few, so larger safe distance
		if (ees[nee_get - 2] <= val)/* this means occupation of last but one ees larger than threshold, so not enough orbitals */ {
			eHbound = ef + (occuMinimumEkBTSecond * kBt > 2 * CyclotronGap) ? occuMinimumEkBTSecond * kBt : 2 * CyclotronGap; // go to the energy upper bounds calculated by current Ef
			return SHIFT_EBOUND_STATUS;// need re-diagonalize
		}
		else {
			int miniFirstNumber = 0;
			for (int it = 0; it < nee_get; ++it) {
				if (ees[it] >= val - valTolerance) { miniFirstNumber = it; break; }
			}
			eHbound = (occuMinimumEkBTSecond - occuMinimumEkBTFirst) * kBt;
			setforwardNumber = (miniFirstNumber - nee_want > 0) ? miniFirstNumber - nee_want : 0;//nee_want may be the Landau level degeneracy, such that can bear with level crossing around EF
			eHbound = ((ees[miniFirstNumber] - ees[setforwardNumber]) > eHbound) ? (ees[miniFirstNumber] - ees[setforwardNumber]) : eHbound;
			//eHbound = (2.0 * CyclotronGap > eHbound) ? 2.0 * CyclotronGap : eHbound;modified from shiftEboudFirst
			eehigh = ees[miniFirstNumber] + eHbound;// use max to make sure the desired miniFirstNumber can still be found next time
			SHIFT_EBOUND_STATUS = true;
			return SHIFT_EBOUND_STATUS;
		}
	}
	else if (nee_get <= 0) {// This is not possible to happen given that we calculate fermiEnergy etc previously
		eHbound = ef + (occuMinimumEkBTSecond * kBt > 2 * CyclotronGap) ? occuMinimumEkBTSecond * kBt : 2 * CyclotronGap; // go to the energy upper bounds calculated by current Ef
		return SHIFT_EBOUND_STATUS;// need re-diagonalize
	}

	return SHIFT_EBOUND_STATUS;
}




/*this (1) examines if the current fermi and occu calculation is reasonable, (2) also dicedes if should re-digonalize.
(3) and make next engenvalue range guess. */
bool shiftEboudFirstOccuTemperatureSorted(int nee_want, int& nee_get, arma::rowvec ees, double& eehigh, double& eelow, double& kBt, double& ef, double* occu, const bool& abool, bool& bbool) {
	bool SHIFT_EBOUND_STATUS = false;
	/*
	three key concepts:
1.	a "band width" defined by the fixed number orbitals contained, e.g., obtained by setbackwardNumber
2.	an energy distance equialence of decay of occupation number, e.g., occuMinimumEkBTFirst * kBt
3.	also the CyclotronGap
	
	for upper bound, apply 2.kBT mostly/even only, prefer not to use 3.CyclotronGap; 
	*/
	
	/*low boundary safe distance can be set by condisering band width from orbital number, temeprature, and energy gap */
	//as as in shiftEboudFirst 
	int setbackwardNumber = (nee_want < nee_get) ? nee_want : nee_get;//pick the small one, usually nee_want is smaller; this is something like a "band width" that contains certain number of orbitals 
	double eLboundNumber = ees[setbackwardNumber - 1] - ees[0];
	double eLboundTemperature = occuMinimumEkBTFirst * kBt;
	double eLbound = (eLboundTemperature > eLboundNumber) ? eLboundTemperature : eLboundNumber;
	eLbound = (CyclotronGap > eLbound) ? CyclotronGap : eLbound;//pick the largest energy distance, usually eLboundNumber is larger?;
	eelow = ees[0] - eLbound;

	/*upper boundary safe distance can be set by condisering orbital number, temeprature, and energy gap, as well as nee_get */
	// prefer to use kBT than in shiftEboudFirst, no CyclotronGap used
	double val = ef + occuMinimumEkBTFirst * kBt;/*This is used to determine an energy value where occupation is eqaul to a threshold Cmin, which exists in my imagination only;
	also is the lower limit for eehigh, additionally higher energies are added to make sure there will be ees higher than val in the next iteration*/
	double valTolerance = 0.001 * kBt;//This is to find ees with energy >= val, but use ees>val-valTolerance in case that energy >= val is missed due to numerical error
	double eHbound;
	int setforwardNumber;//This is also used to calculate something like a "band width" defined by the fixed number orbitals contained
	/*Notice we can allow recalculating occupationNumber and fermiEnergy here by using only a desired number of orbitals, this is useful because occu may very sensitively depend on orbital number and
	I do not want to calculate that many of them, so I fix it.*/
	int maxNumOrbital = 1000;
	if (nee_get > maxNumOrbital) {// here is a situation where the orbital number is too many, so smaller safe distance
		if (ees[maxNumOrbital - 1] <= val)/* this means occupation of ees[999=maxNumOrbital-1] larger than threshold Cmin */ {
			//eHbound = (occuMinimumEkBTSecond - occuMinimumEkBTFirst) * kBt; modified from shiftEboudFirst
			setforwardNumber = (maxNumOrbital - 100 > 0) ? maxNumOrbital - 100 : 0;
			eHbound = ees[maxNumOrbital - 1] - ees[setforwardNumber];//make sure have ~100 more than maxNumOrbital
			//eHbound = ((ees[maxNumOrbital - 1] - ees[setforwardNumber]) > eHbound) ? (ees[maxNumOrbital - 1] - ees[setforwardNumber]) : eHbound;modified from shiftEboudFirst
			//eHbound = (0.5 * CyclotronGap > eHbound) ? 0.5 * CyclotronGap : eHbound; modified from shiftEboudFirst
			if (abool) {
				double deltaNdiff; int icount = -1;//this definition is only to make up for the reuiqred number of arguments for findOccupationNumber
				if (findOccupationNumberSorted(nee_want, maxNumOrbital, ees, eehigh, eelow, kBt, ef, occu, deltaNdiff, icount)) {//This means the orbitals more than maxNumOrbital are ignored
					eehigh = ees[maxNumOrbital - 1] + eHbound;
					nee_get = maxNumOrbital;// use max to make sure smooth evolution of the relation nee_get > maxNumOrbital, such that do not miss/lose orbitals suddenly next time 
					bbool = true;
					SHIFT_EBOUND_STATUS = true;
					return SHIFT_EBOUND_STATUS;
					/*PHILOSOPHY: new ef_new will be larger than old ef_old due to less orbitals, so if val_old=ef_old+xx is larger than es[maxNumOrbital - 1], 
					then val_new=ef_new+xx >= es[maxNumOrbital - 1];
					*in next solutions, orbital numbers be ~maxNumOrbital+100, ef_next_iter_old be around last ef_new or less since ~maxNumOrbital+100> maxNumOrbital, or
					* val_old < ef_next_iter_old < val_new, this would risk ef_next_iter_old+xx >= es[maxNumOrbital - 1] and 
					as a result go to the following if-loop branch: ef_next_iter_old+xx >= es[maxNumOrbital - 1];
					which will then involve kBT to set upper limit to include even larger # orbitals beyond ~maxNumOrbital+100; 
					To walk around this problem: use the smaller upper limit that are determined from kBT and setforwardNumber in the other if-loop branch
					REASON: in the next if-loop branch, ef_next_iter_old+xx is larger than it "should be when orbital # infinite", thus eHbound_kBT = ef + occuMinimumEkBTSecond * kBt is overestimated,
					if eHbound_kBT < eHbound_setforwardNumber, then it means ees[~maxNumOrbital+100] > Cmin when the orbital # constraint of ~maxNumOrbital+100 is released.
					*/
				}
				else {
					eHbound = (1.0 * CyclotronGap > eHbound) ? 1.0 * CyclotronGap : eHbound;//modified from shiftEboudFirst
					eehigh = ees[nee_get - 1] + eHbound;// This choice happens only if above does not find fermi energy, does not seem possible; Expand search region in this case;
					return SHIFT_EBOUND_STATUS; // need re-diagonalize
				}
			}
		}
		else /* this means occupation of ees[999=maxNumOrbital-1] smaller than threshold Cmin */ {//"completely" modified from shiftEboudFirst
			eHbound = ef + occuMinimumEkBTSecond * kBt - ees[maxNumOrbital - 1];//this might be negative
			setforwardNumber = (maxNumOrbital-1 - 100 > 0) ? maxNumOrbital-1 - 100 : 0;
			eHbound = ((ees[maxNumOrbital - 1] - ees[setforwardNumber]) < eHbound) ? (ees[maxNumOrbital - 1] - ees[setforwardNumber]) : eHbound;//pick the smaller one
			//eHbound = (0.5 * CyclotronGap > eHbound) ? 0.5 * CyclotronGap : eHbound; modified from shiftEboudFirst
			eehigh = ees[maxNumOrbital - 1] + eHbound;// use min to make sure the orbital # can decrease when it "should", and can not increase when it "should not"
			SHIFT_EBOUND_STATUS = true;
			return SHIFT_EBOUND_STATUS;
		}
	}
	else if (nee_get > 3 * nee_want && nee_get <= maxNumOrbital) {// here is a situation where the orbital number is not too many but not too few, , so moderate safe distance
		if (ees[nee_get - 2] <= val)/* this means occupation of last but one ees larger than threshold Cmin, so not enough orbitals */ {
			eHbound = ef + (occuMinimumEkBTSecond * kBt > 1.0 * CyclotronGap) ? occuMinimumEkBTSecond * kBt : 1.0 * CyclotronGap; // go to the energy upper bounds calculated by current Ef
			return SHIFT_EBOUND_STATUS;// need re-diagonalize
		}
		else {
			int miniFirstNumber = 0;
			for (int it = 0; it < nee_get; ++it) {
				if (ees[it] >= val - valTolerance) { miniFirstNumber = it; break; }//find the first/lowest energy orbital to have occupation less than Cmin
			}
			eHbound = (occuMinimumEkBTSecond - occuMinimumEkBTFirst) * kBt;
			setforwardNumber = (miniFirstNumber - nee_want > 0) ? miniFirstNumber - nee_want : 0;//nee_want may be the Landau level degeneracy, such that can bear with level crossing around EF
			eHbound = ((ees[miniFirstNumber] - ees[setforwardNumber]) > eHbound) ? (ees[miniFirstNumber] - ees[setforwardNumber]) : eHbound;
			//eHbound = (0.5 * CyclotronGap > eHbound) ? 0.5 * CyclotronGap : eHbound; modified from shiftEboudFirst
			eehigh = ees[miniFirstNumber] + eHbound;// use max to make sure the desired miniFirstNumber can still be found next time
			SHIFT_EBOUND_STATUS = true;
			return SHIFT_EBOUND_STATUS;
		}
	}
	else if (nee_get <= 3 * nee_want && nee_get >= 1) {// here is a situation where the orbital number is too few, so larger safe distance
		if (ees[nee_get - 2] <= val)/* this means occupation of last but one ees larger than threshold, so not enough orbitals */ {
			eHbound = ef + (occuMinimumEkBTSecond * kBt > 2 * CyclotronGap) ? occuMinimumEkBTSecond * kBt : 2 * CyclotronGap; // go to the energy upper bounds calculated by current Ef
			return SHIFT_EBOUND_STATUS;// need re-diagonalize
		}
		else {
			int miniFirstNumber = 0;
			for (int it = 0; it < nee_get; ++it) {
				if (ees[it] >= val - valTolerance) { miniFirstNumber = it; break; }
			}
			eHbound = (occuMinimumEkBTSecond - occuMinimumEkBTFirst) * kBt;
			setforwardNumber = (miniFirstNumber - nee_want > 0) ? miniFirstNumber - nee_want : 0;//nee_want may be the Landau level degeneracy, such that can bear with level crossing around EF
			eHbound = ((ees[miniFirstNumber] - ees[setforwardNumber]) > eHbound) ? (ees[miniFirstNumber] - ees[setforwardNumber]) : eHbound;
			//eHbound = (2.0 * CyclotronGap > eHbound) ? 2.0 * CyclotronGap : eHbound;modified from shiftEboudFirst
			eehigh = ees[miniFirstNumber] + eHbound;// use max to make sure the desired miniFirstNumber can still be found next time
			SHIFT_EBOUND_STATUS = true;
			return SHIFT_EBOUND_STATUS;
		}
	}
	else if (nee_get <= 0) {// This is not possible to happen given that we calculate fermiEnergy etc previously
		eHbound = ef + (occuMinimumEkBTSecond * kBt > 2 * CyclotronGap) ? occuMinimumEkBTSecond * kBt : 2 * CyclotronGap; // go to the energy upper bounds calculated by current Ef
		return SHIFT_EBOUND_STATUS;// need re-diagonalize
	}

	return SHIFT_EBOUND_STATUS;
}




bool shiftEboudFirstOccuTemperatureOnly(int nee_want, int& nee_get, const double* ees, double& eehigh, double& eelow, double& kBt, double& ef, double* occu, const bool& abool, bool& bbool) {
	bool SHIFT_EBOUND_STATUS = false;
	/*
	three key concepts:
1.	a "band width" defined by the fixed number orbitals contained, e.g., obtained by setbackwardNumber
2.	an energy distance equialence of decay of occupation number, e.g., occuMinimumEkBTFirst * kBt
3.	also the CyclotronGap

	for upper bound, apply 2.kBT mostly/even only, prefer not to use 3.CyclotronGap
	*/

	/*low boundary safe distance can be set by condisering band width from orbital number, temeprature, and energy gap */
	//as as in shiftEboudFirst 
	int setbackwardNumber = (nee_want < nee_get) ? nee_want : nee_get;//pick the small one, usually nee_want is smaller; this is something like a "band width" that contains certain number of orbitals 
	double eLboundNumber = ees[setbackwardNumber - 1] - ees[0];
	double eLboundTemperature = occuMinimumEkBTFirst * kBt;
	double eLbound = (eLboundTemperature > eLboundNumber) ? eLboundTemperature : eLboundNumber;
	eLbound = (CyclotronGap > eLbound) ? CyclotronGap : eLbound;//pick the largest energy distance, usually eLboundNumber is larger?;
	eelow = ees[0] - eLbound;

	/*upper boundary safe distance can be set by condisering orbital number, temeprature, and energy gap, as well as nee_get */
	// prefer to use kBT than in shiftEboudFirst, no CyclotronGap used
	double val = ef + occuMinimumEkBTFirst * kBt;/*This is used to determine an energy value where occupation is eqaul to a threshold Cmin, which exists in my imagination only;
	also is the lower limit for eehigh, additionally higher energies are added to make sure there will be ees higher than val in the next iteration*/
	double valTolerance = 0.001 * kBt;//This is to find ees with energy >= val, but use ees>val-valTolerance in case that energy >= val is missed due to numerical error
	double eHbound;
	int setforwardNumber;//This is also used to calculate something like a "band width" defined by the fixed number orbitals contained
	/*Notice we can allow recalculating occupationNumber and fermiEnergy here by using only a desired number of orbitals, this is useful because occu may very sensitively depend on orbital number and
	I do not want to calculate that many of them, so I fix it.*/
	int maxNumOrbital = 1000;
	if (nee_get > maxNumOrbital) {// here is a situation where the orbital number is too many, so smaller safe distance
		if (ees[maxNumOrbital - 1] <= val)/* this means occupation of ees[999=maxNumOrbital-1] larger than threshold Cmin */ {
			//eHbound = (occuMinimumEkBTSecond - occuMinimumEkBTFirst) * kBt; modified from shiftEboudFirst
			setforwardNumber = (maxNumOrbital - 100 > 0) ? maxNumOrbital - 100 : 0;
			eHbound = ees[maxNumOrbital - 1] - ees[setforwardNumber];//make sure have ~100 more than maxNumOrbital
			//eHbound = ((ees[maxNumOrbital - 1] - ees[setforwardNumber]) > eHbound) ? (ees[maxNumOrbital - 1] - ees[setforwardNumber]) : eHbound;modified from shiftEboudFirst
			//eHbound = (0.5 * CyclotronGap > eHbound) ? 0.5 * CyclotronGap : eHbound; modified from shiftEboudFirst
			if (abool) {
				double deltaNdiff; int icount = -1;//this definition is only to make up for the reuiqred number of arguments for findOccupationNumber
				if (findOccupationNumber(nee_want, maxNumOrbital, ees, eehigh, eelow, kBt, ef, occu, deltaNdiff, icount)) {//This means the orbitals more than maxNumOrbital are ignored
					eehigh = ees[maxNumOrbital - 1] + eHbound;
					nee_get = maxNumOrbital;// use max to make sure smooth evolution of the relation nee_get > maxNumOrbital, such that do not miss/lose orbitals suddenly next time 
					bbool = true;
					SHIFT_EBOUND_STATUS = true;
					return SHIFT_EBOUND_STATUS;
					/*PHILOSOPHY: new ef_new will be larger than old ef_old due to less orbitals, so if val_old=ef_old+xx is larger than es[maxNumOrbital - 1],
					then val_new=ef_new+xx >= es[maxNumOrbital - 1];
					*in next solutions, orbital numbers be ~maxNumOrbital+100, ef_next_iter_old be around last ef_new or less since ~maxNumOrbital+100> maxNumOrbital, or
					* val_old < ef_next_iter_old < val_new, this would risk ef_next_iter_old+xx >= es[maxNumOrbital - 1] and
					as a result go to the following if-loop branch: ef_next_iter_old+xx >= es[maxNumOrbital - 1];
					which will then involve kBT to set upper limit to include even larger # orbitals beyond ~maxNumOrbital+100;
					To walk around this problem: use the smaller upper limit that are determined from kBT and setforwardNumber in the other if-loop branch
					REASON: in the next if-loop branch, ef_next_iter_old+xx is larger than it "should be when orbital # infinite", thus eHbound_kBT = ef + occuMinimumEkBTSecond * kBt is overestimated,
					if eHbound_kBT < eHbound_setforwardNumber, then it means ees[~maxNumOrbital+100] > Cmin when the orbital # constraint of ~maxNumOrbital+100 is released.
					*/
				}
				else {
					eHbound = (1.0 * CyclotronGap > eHbound) ? 1.0 * CyclotronGap : eHbound;//modified from shiftEboudFirst
					eehigh = ees[nee_get - 1] + eHbound;// This choice happens only if above does not find fermi energy, does not seem possible; Expand search region in this case;
					return SHIFT_EBOUND_STATUS; // need re-diagonalize
				}
			}
		}
		else /* this means occupation of ees[999=maxNumOrbital-1] smaller than threshold Cmin */ {//"completely" modified from shiftEboudFirst
			eHbound = ef + occuMinimumEkBTSecond * kBt - ees[maxNumOrbital - 1];//this might be negative
			setforwardNumber = (maxNumOrbital - 1 - 100 > 0) ? maxNumOrbital - 1 - 100 : 0;
			eHbound = ((ees[maxNumOrbital - 1] - ees[setforwardNumber]) < eHbound) ? (ees[maxNumOrbital - 1] - ees[setforwardNumber]) : eHbound;//pick the smaller one
			//eHbound = (0.5 * CyclotronGap > eHbound) ? 0.5 * CyclotronGap : eHbound; modified from shiftEboudFirst
			eehigh = ees[maxNumOrbital - 1] + eHbound;// use min to make sure the orbital # can decrease when it "should", and can not increase when it "should not"
			SHIFT_EBOUND_STATUS = true;
			return SHIFT_EBOUND_STATUS;
		}
	}
	else if (nee_get > 3 * nee_want && nee_get <= maxNumOrbital) {// here is a situation where the orbital number is not too many but not too few, , so moderate safe distance
		if (ees[nee_get - 2] <= val)/* this means occupation of last but one ees larger than threshold Cmin, so not enough orbitals */ {
			eHbound = ef + (occuMinimumEkBTSecond * kBt > 1.0 * CyclotronGap) ? occuMinimumEkBTSecond * kBt : 1.0 * CyclotronGap; // go to the energy upper bounds calculated by current Ef
			return SHIFT_EBOUND_STATUS;// need re-diagonalize
		}
		else {
			int miniFirstNumber = 0;
			for (int it = 0; it < nee_get; ++it) {
				if (ees[it] >= val - valTolerance) { miniFirstNumber = it; break; }//find the first/lowest energy orbital to have occupation less than Cmin
			}
			eHbound = (occuMinimumEkBTSecond - occuMinimumEkBTFirst) * kBt;
			setforwardNumber = (miniFirstNumber - 50 > 0) ? miniFirstNumber - 50 : 0;
			eHbound = ((ees[miniFirstNumber] - ees[setforwardNumber]) > eHbound) ? (ees[miniFirstNumber] - ees[setforwardNumber]) : eHbound;
			//eHbound = (0.5 * CyclotronGap > eHbound) ? 0.5 * CyclotronGap : eHbound; modified from shiftEboudFirst
			eehigh = ees[miniFirstNumber] + eHbound;// use max to make sure the desired miniFirstNumber can still be found next time
			SHIFT_EBOUND_STATUS = true;
			return SHIFT_EBOUND_STATUS;
		}
	}
	else if (nee_get <= 3 * nee_want && nee_get >= 1) {// here is a situation where the orbital number is too few, so larger safe distance
		if (ees[nee_get - 2] <= val)/* this means occupation of last but one ees larger than threshold, so not enough orbitals */ {
			eHbound = ef + (occuMinimumEkBTSecond * kBt > 2 * CyclotronGap) ? occuMinimumEkBTSecond * kBt : 2 * CyclotronGap; // go to the energy upper bounds calculated by current Ef
			return SHIFT_EBOUND_STATUS;// need re-diagonalize
		}
		else {
			int miniFirstNumber = 0;
			for (int it = 0; it < nee_get; ++it) {
				if (ees[it] >= val - valTolerance) { miniFirstNumber = it; break; }
			}
			eHbound = (occuMinimumEkBTSecond - occuMinimumEkBTFirst) * kBt;
			setforwardNumber = (miniFirstNumber - 50 > 0) ? miniFirstNumber - 50 : 0;
			eHbound = ((ees[miniFirstNumber] - ees[setforwardNumber]) > eHbound) ? (ees[miniFirstNumber] - ees[setforwardNumber]) : eHbound;
			//eHbound = (2.0 * CyclotronGap > eHbound) ? 2.0 * CyclotronGap : eHbound;modified from shiftEboudFirst
			eehigh = ees[miniFirstNumber] + eHbound;// use max to make sure the desired miniFirstNumber can still be found next time
			SHIFT_EBOUND_STATUS = true;
			return SHIFT_EBOUND_STATUS;
		}
	}
	else if (nee_get <= 0) {// This is not possible to happen given that we calculate fermiEnergy etc previously
		eHbound = ef + (occuMinimumEkBTSecond * kBt > 2 * CyclotronGap) ? occuMinimumEkBTSecond * kBt : 2 * CyclotronGap; // go to the energy upper bounds calculated by current Ef
		return SHIFT_EBOUND_STATUS;// need re-diagonalize
	}

	return SHIFT_EBOUND_STATUS;
}


//calibrateEigenBoundsTemperature(Np, Nee, Eigenvals, Energyhigh, Energylow, kBT, fermiEnergy, occupationNumber, Nsubspace, estimateDos) != 0
//yayun 240130
int calibrateFixEigenSorted(int nee_want, int& nee_get, arma::rowvec ees, double& eehigh, double& eelow, double& kBt, double& Ef, double* occu, int& nee_expect, 
	double& averDos, 
	int& nee_fix, int fmode, const int& nrun, const double& PotentialMin, double& lastemin, double& nee_wantdecimal) {
	using namespace std;
	using namespace MAIN;
	const string	filefermiNotFound = prefix + "NotFound.sum";//This file is uesd to contain & flag unexpected errors when determining fermi energy
	int ret = 0;
	double last;
	double OldEnergyWindowWidth = eehigh -eelow;
	averDos = double(nee_get + 5.0) / (ees[nee_get - 1] - ees[0]); // maximum estimate of average density of states using output energy

	if (fmode==0) {//use all orbitals
		//find occupation Number and fermi energy
		double deltaNdiff; int icount = -1;// this is the numerical difference between total occupied and the integer nee_want, should be small
		if (findOccupationNumberSorted(nee_want, nee_get, ees, eehigh, eelow, kBt, Ef, occu, deltaNdiff, icount)) {
			//reportOccupation(-1, nee_get);
			mainout << "This is the temperature test, fermiEnergy=" << Ef << "\t" << "Old Energy bound guesses: " << std::setprecision(6) << eelow << " -- " << eehigh << endl;
		}
		else {
			mainout << "Fail to find fermi energy, fermiEnergy=" << Ef << ",icount is " << icount << ",deltaNdiff is" << std::setprecision(18) << deltaNdiff << "\t" << "Old Energy bound guesses: " << std::setprecision(6) << eelow << " -- " << eehigh << endl;
			eelow = ees[0] - CyclotronGap;
			eehigh = ees[nee_get - 1] + CyclotronGap;
			//nee_expect = min(int((eehigh - eelow) * averDos * 1.5), MaxSubspace - 10);//expect new Norb using existing dos
			nee_expect = min(int((eehigh - eelow - OldEnergyWindowWidth) * averDos * 1.5) + nee_get, MaxSubspace - 10);//expect new Norb using existing dos and existing orbitals
			ret = 1;
			return ret;
		}
		//change the bounds. Also possibly change occupationNumber and fermiEnergy!!!!
		bool modifyFermiOccu = true;// This is to say when there are too many orbitals, then use only certain number of them to calculate ef and occu
		bool modified = false;//tell mainout if ef and occu modified
		if (shiftEboudFirstOccuTemperatureSorted(nee_want, nee_get, ees, eehigh, eelow, kBt, Ef, occu, modifyFermiOccu, modified)) {
			mainout << "New energy bounds set=" << "\t" << "New Energy bound guesses: " << eelow << " -- " << eehigh << endl;
			if (modified) {
				//reportOccupation(-2, nee_get);
				mainout << "New occupations are reset as well as fermiEnergy, which is=" << Ef << endl;
			}
			//nee_expect = min(int((eehigh - eelow) * averDos * 1.5), MaxSubspace - 10);
			nee_expect = min(int((eehigh - eelow - OldEnergyWindowWidth) * averDos * 1.5) + nee_get, MaxSubspace - 10);//expect new Norb using existing dos and existing orbitals
			/*this line below is added on 20210209, to use all orbitals from diagonalization, to be consistent with nee_fix*/
			nee_fix = nee_get;
			return ret;
		}
		else {
			mainout << "Fail to change bounds, need new bounds to re-diagonalize =" << Ef << "\t" << "New Energy bound guesses: " << eelow << " -- " << eehigh << endl;
			//nee_expect = min(int((eehigh - eelow) * averDos * 1.5), MaxSubspace - 10);
			nee_expect = min(int((eehigh - eelow - OldEnergyWindowWidth) * averDos * 1.5) + nee_get, MaxSubspace - 10);//expect new Norb using existing dos and existing orbitals
			ret = 1;
			return ret;
		}
	}else if(fmode == 1) {//use up to the first orbital that is occupied less than C_min (see syspara.hpp for def); here operations can be crude, real fix later;
		//find occupation Number and fermi energy
		double deltaNdiff; // this is the numerical difference between total occupied and the integer nee_want, should be small
		int icount = -1;//is the number of search in the fermi - distribution function to find the fermi energy		
		lastemin = ees[0]; nee_wantdecimal = nee_want;//here the particle number is fixed, so these values are known and given like this, these data not used, just for storage;
		//int FERMI_ENERGY_STATUS_NUM = findOccuCminLocateOneaboveEf(nee_want, nee_get, ees, eehigh, eelow, kBt, Ef, occu, deltaNdiff, icount, nee_fix, PotentialMin);
		
		int nee_expect_temp;
		int FERMI_ENERGY_STATUS_NUM = findOccukBtGuessRange_Sorted(nee_expect_temp, nee_want, nee_get, ees, eehigh, eelow, kBt, Ef, occu, deltaNdiff, icount, nee_fix, PotentialMin);
		if (FERMI_ENERGY_STATUS_NUM == 3 || FERMI_ENERGY_STATUS_NUM == 4) {
			nee_expect = min(int((eehigh - eelow - OldEnergyWindowWidth) * averDos * 1.5) + nee_get, MaxSubspace - 10);//expect new Norb using existing dos and existing orbitals
			nee_expect = max(nee_expect_temp, nee_expect);//yayun, 20240127, it happens that Nsubspace is too small when FERMI_ENERGY_STATUS_NUM == 4 "E:\sharedata01\yhu\workspace\data\CFDFTsquare\LXfracFixNee240116\28\slurm-35187_1.out"
			mainout << "Found: " << FERMI_ENERGY_STATUS_NUM << " fermiEnergy= " << Ef << "\t" << "New bound guesses: " << std::setprecision(6) << eelow << " -- " << eehigh << endl;									
			return ret;//next iteration
		}
		else if (FERMI_ENERGY_STATUS_NUM == 5) {
			//nee_expect = min(int((eehigh - eelow) * averDos * 1.5), MaxSubspace - 10);//expect new Norb using existing dos
			nee_expect = min(int((eehigh - eelow - OldEnergyWindowWidth) * averDos * 1.5) + nee_get, MaxSubspace - 10);//expect new Norb using existing dos and existing orbitals
			mainout << "Special, large fermi when nee_want=nee_get, Found: FERMI_ENERGY_STATUS_NUM = " << FERMI_ENERGY_STATUS_NUM << " fermiEnergy=" << Ef << "\t" << "New bound guesses: " << std::setprecision(6) << eelow << " -- " << eehigh << endl;
			return ret;//next iteration
		}
		else if(FERMI_ENERGY_STATUS_NUM == 1){
			mainout << "Fail to find Ef, FERMI_ENERGY_STATUS_NUM = "<< FERMI_ENERGY_STATUS_NUM <<", fermiEnergy= " << Ef << ", icount is " << icount << ", deltaNdiff is " << std::setprecision(18) << deltaNdiff << "\t" << "New bound guesses: " << std::setprecision(6) << eelow << " -- " << eehigh << endl;
			//nee_expect = min(int((eehigh - eelow) * averDos * 1.5), MaxSubspace - 10);//expect new Norb using existing dos	
			nee_expect = min(int((eehigh - eelow - OldEnergyWindowWidth) * averDos * 1.5) + nee_get, MaxSubspace - 10);//expect new Norb using existing dos and existing orbitals
			//nrun,nee_get,nee_expect,nee_fix,icount,deltaNdiff,eehigh,eelow,Ef
			fermiNotFound(filefermiNotFound,nrun,nee_get,nee_expect,nee_fix,icount,deltaNdiff,eehigh,eelow,Ef,FERMI_ENERGY_STATUS_NUM);
			ret = 1;
			return ret;//rediagonalize for fermi energy
		}
		else if (FERMI_ENERGY_STATUS_NUM == 2) {
			//nee_expect = min(int((eehigh - eelow) * averDos * 1.5), MaxSubspace - 10);//expect new Norb using existing dos
			nee_expect = min(int((eehigh - eelow - OldEnergyWindowWidth) * averDos * 1.5) + nee_get, MaxSubspace - 10);//expect new Norb using existing dos and existing orbitals
			mainout << "Cmin not satisfied, FERMI_ENERGY_STATUS_NUM= " << FERMI_ENERGY_STATUS_NUM << " eigen low&high: " << ees[0] <<"--" << ees[nee_want - 1] <<  ",fermiEnergy= " << Ef << " New bound guesses: " << std::setprecision(6) << eelow << " -- " << eehigh << endl;
			fermiNotFound(filefermiNotFound, nrun, nee_get, nee_expect, nee_fix, icount, deltaNdiff, eehigh, eelow, Ef, FERMI_ENERGY_STATUS_NUM);
			ret = 1;
			return ret;//rediagonalize for a few more orbitals
		}		
	}
	else if (fmode == 2) {//find occupation Number
		int FERMI_ENERGY_STATUS_NUM = findOccuCminfixedEfSorted(nee_want, nee_get, ees, eehigh, eelow, kBt, Ef, occu, nee_fix, PotentialMin, lastemin, nee_wantdecimal);
		if (FERMI_ENERGY_STATUS_NUM == 2) {//this condition is the only possibility
			//nee_expect = min(int((eehigh - eelow) * averDos * 1.5), MaxSubspace - 10);//expect new Norb using existing dos
			nee_expect = min(int((eehigh - eelow - OldEnergyWindowWidth) * averDos * 1.5) + nee_get, MaxSubspace - 10);//expect new Norb using existing dos and existing orbitals
			mainout << std::setprecision(6) << " Emin: " << lastemin << "\t" << "New bound guesses: " << eelow << " -- " << eehigh << endl;
			return ret;//next iteration
		}
	}
	return ret;//will not happen
}




/*
 * Simple checks to see whether enough eigenstates are found and sets new eigenbounds.
 * Also checks if lower bound is low enough (Need nee_get > 1), assuming E is bounded below.
 * If eigensolver errs, in which case nee_get == -1, noop and return -1.
 * Return 0 if the checks are passed. Nonzero otherwise.
 * In both cases the eigenbounds are set.
 */
int calibrateEigenBounds(int nee_want, int nee_get, const double *ees, double &eehigh, double &eelow) {
	using namespace std;
	using namespace MAIN;
	int ret = 0;
	double last;

	mainout << "Current energy bounds:\t"<<Energylow<<" -- "<<Energyhigh<<endl;
	mainout << "Current eigenval range:\t"<<ees[0]<<" -- ";
	if(nee_want <= nee_get) {
		mainout << ees[nee_want-1] << endl;
	} else {
		mainout << "NaN (" << ees[nee_get] << ")" << endl;
	}

	
	if(nee_want < 2) {
		auto msg = "WARN: calibrateEigenBounds: less than 2 states dersired.";
		mainout << msg << endl;
	}

	/* check if lower bound is low enough */
	if (nee_get > 1) {
		if(ees[0]-eelow < ees[1]-ees[0]) {
			ret += 1;
			last = eelow;
			eelow -= 10*(ees[1] - ees[0]);
			mainout << "Need to decrease Eigenlow, old: "<<last<<", new:"<<eelow<<endl;
		} else if(ees[0]-eelow > 3*(ees[1]-ees[0])) {
			last = eelow;
			eelow = ees[0] - max(3*(ees[1]-ees[0]),0.5*CyclotronGap);
			mainout << "Will increase Eigenlow, old: "<<last<<", new:"<<eelow<<endl;
		}
	}

	/*
	 * Check & adjust upper bound
	 * It is assume that during recusion, nee_get < nee_want occurs mainly due to eehigh being too low.
	 */
	if(nee_get < nee_want) {
		last = eehigh;
		if(nee_get > 5) {
			eehigh+= ((ees[nee_get-1] - ees[nee_get-5]) / 5.0) * (nee_want-nee_get)*1.3;
		} else {
			eehigh+= CyclotronGap * (nee_want-nee_get) / 1.0 / MAIN::Nstate_per_Elevel;
		}
		ret += 1;
		mainout << "Need to increase Eigenhigh, old: "<<last<<", new:"<<eehigh<<endl;
	} else if(nee_get > nee_want * 1.5) {
		last = eehigh;
		eehigh+= ((ees[nee_get-1] - ees[nee_get-5]) / 5.0) * (nee_want-nee_get)/2.0;
		mainout << "Will decrease Eigenhigh, old: "<<last<<", new:"<<eehigh<<endl;
	}

	return ret;
}



int calibrateEigenBoundsTemperature(int nee_want, int& nee_get, const double* ees, double& eehigh, double& eelow, double& kBt, double& Ef, double *occu, int& nee_expect, double& averDos) {
	using namespace std;
	using namespace MAIN;
	int ret = 0;
	double last;
	averDos = double(nee_get+5.0)/(ees[nee_get - 1] - ees[0]); // maximum estimate of density of states using output energy

	if (isThereTemperature == false) {

		mainout << "Current energy bounds:\t" << Energylow << " -- " << Energyhigh << endl;
		mainout << "Current eigenval range:\t" << ees[0] << " -- ";
		if (nee_want <= nee_get) {
			mainout << ees[nee_want - 1] << endl;
		}

		else {
			mainout << "NaN (" << ees[nee_get] << ")" << endl;
		}


		if (nee_want < 2) {
			auto msg = "WARN: calibrateEigenBounds: less than 2 states dersired.";
			mainout << msg << endl;
		}

		/* check if lower bound is low enough */
		if (nee_get > 1) {
			if (ees[0] - eelow < ees[1] - ees[0]) {
				ret += 1;
				last = eelow;
				eelow -= 10 * (ees[1] - ees[0]);
				mainout << "Need to decrease Eigenlow, old: " << last << ", new:" << eelow << endl;
			}
			else if (ees[0] - eelow > 3 * (ees[1] - ees[0])) {
				last = eelow;
				eelow = ees[0] - max(3 * (ees[1] - ees[0]), 0.5 * CyclotronGap);
				mainout << "Will increase Eigenlow, old: " << last << ", new:" << eelow << endl;
			}
		}

		/*
		 * Check & adjust upper bound
		 * It is assume that during recusion, nee_get < nee_want occurs mainly due to eehigh being too low.
		 */
		if (nee_get < nee_want) {
			last = eehigh;
			if (nee_get > 5) {
				eehigh += ((ees[nee_get - 1] - ees[nee_get - 5]) / 5.0) * (nee_want - nee_get) * 1.3;
			}
			else {
				eehigh += CyclotronGap * (nee_want - nee_get) / 1.0 / MAIN::Nstate_per_Elevel;
			}
			ret += 1;
			mainout << "Need to increase Eigenhigh, old: " << last << ", new:" << eehigh << endl;
		}
		else if (nee_get > nee_want * 1.5) {
			last = eehigh;
			eehigh += ((ees[nee_get - 1] - ees[nee_get - 5]) / 5.0) * (nee_want - nee_get) / 2.0;
			mainout << "Will decrease Eigenhigh, old: " << last << ", new:" << eehigh << endl;
		}		
		return ret;
		
	}else{		
		//find occupation Number and fermi energy
		double deltaNdiff = -1.0; int icount = -1;// this is the numerical difference between total occupied and the integer nee_want, should be small
		if (findOccupationNumber(nee_want, nee_get, ees, eehigh, eelow, kBt, Ef, occu, deltaNdiff, icount)) {
			//reportOccupation(-1, nee_get);
			mainout << "This is the temperature test, fermiEnergy=" << Ef << "\t" << "Old Energy bound guesses: " << std::setprecision(6) << eelow << " -- " << eehigh << endl;
		}
		else {
			//mainout << "deltaN:" << deltaNdiff << ", icount:" << icount << endl;
			mainout << "Fail to find Ef, fermiEnergy=" << Ef << ",icount is " << icount << ",deltaNdiff is " << std::setprecision(18) << deltaNdiff << "\t" << "Old Energy bound guesses: " << std::setprecision(6) << eelow << " -- " << eehigh << endl;
			eelow = ees[0] - CyclotronGap;
			eehigh = ees[nee_get - 1] + CyclotronGap;
			nee_expect = min(int((eehigh - eelow) * averDos*1.5),MaxSubspace-10);//expect new Norb using existing dos
			ret = 1;
			return ret;
		}
		//change the bounds. Also possibly change occupationNumber and fermiEnergy!!!!
		bool modifyFermiOccu = true;// This is to say when there are too many orbitals, then use only certain number of them to calculate ef and occu
		bool modified = false;//tell mainout if ef and occu modified
		if (shiftEboudFirstOccuTemperature(nee_want, nee_get, ees, eehigh, eelow, kBt, Ef, occu, modifyFermiOccu, modified)) {//shiftEboudFirstOccuTemperature
			mainout << "New energy bounds set=" << "\t" << "New Energy bound guesses: " << eelow << " -- " << eehigh << endl;
			if (modified) {
				//reportOccupation(-2, nee_get);
				mainout << "New occupations are reset as well as fermiEnergy, which is=" << Ef << endl;
			}
			nee_expect = min(int((eehigh - eelow) * averDos * 1.5), MaxSubspace - 10);
			return ret;
		}
		else {
			mainout << "Fail to change bounds, need new bounds to re-diagonalize =" << Ef << "\t" << "New Energy bound guesses: " << eelow << " -- " << eehigh << endl;
			nee_expect = min(int((eehigh - eelow) * averDos * 1.5), MaxSubspace - 10);
			ret = 1;
			return ret;
		}
	}
}






//calibrateEigenBoundsTemperature(Np, Nee, Eigenvals, Energyhigh, Energylow, kBT, fermiEnergy, occupationNumber, Nsubspace, estimateDos) != 0

int calibrateFixEigenZeroth(int nee_want, int& nee_get, const double* ees, double& eehigh, double& eelow, double& kBt, double& Ef, double* occu, int& nee_expect, 
	double& averDos, 
	int& nee_fix, int fmode, const int& nrun, const double& PotentialMin, double& lastemin, double& nee_wantdecimal) {
	using namespace std;
	using namespace MAIN;
	const string	filefermiNotFound = prefix + "NotFound.sum";//This file is uesd to contain & flag unexpected errors when determining fermi energy
	int ret = 0;
	double last;
	averDos = double(nee_get + 5.0) / (ees[nee_get - 1] - ees[0]); // maximum estimate of average density of states using output energy

	if (fmode==0) {//use all orbitals
		//find occupation Number and fermi energy
		double deltaNdiff; int icount = -1;// this is the numerical difference between total occupied and the integer nee_want, should be small
		if (findOccupationNumber(nee_want, nee_get, ees, eehigh, eelow, kBt, Ef, occu, deltaNdiff, icount)) {
			//reportOccupation(-1, nee_get);
			mainout << "This is the temperature test, fermiEnergy=" << Ef << "\t" << "Old Energy bound guesses: " << std::setprecision(6) << eelow << " -- " << eehigh << endl;
		}
		else {
			mainout << "Fail to find fermi energy, fermiEnergy=" << Ef << ",icount is " << icount << ",deltaNdiff is" << std::setprecision(18) << deltaNdiff << "\t" << "Old Energy bound guesses: " << std::setprecision(6) << eelow << " -- " << eehigh << endl;
			eelow = ees[0] - CyclotronGap;
			eehigh = ees[nee_get - 1] + CyclotronGap;
			nee_expect = min(int((eehigh - eelow) * averDos * 1.5), MaxSubspace - 10);//expect new Norb using existing dos
			ret = 1;
			return ret;
		}
		//change the bounds. Also possibly change occupationNumber and fermiEnergy!!!!
		bool modifyFermiOccu = true;// This is to say when there are too many orbitals, then use only certain number of them to calculate ef and occu
		bool modified = false;//tell mainout if ef and occu modified
		if (shiftEboudFirstOccuTemperature(nee_want, nee_get, ees, eehigh, eelow, kBt, Ef, occu, modifyFermiOccu, modified)) {//shiftEboudFirstOccuTemperature
			mainout << "New energy bounds set=" << "\t" << "New Energy bound guesses: " << eelow << " -- " << eehigh << endl;
			if (modified) {
				//reportOccupation(-2, nee_get);
				mainout << "New occupations are reset as well as fermiEnergy, which is=" << Ef << endl;
			}
			nee_expect = min(int((eehigh - eelow) * averDos * 1.5), MaxSubspace - 10);
			/*this line below is added on 20210209, to use all orbitals from diagonalization, to be consistent with nee_fix*/
			nee_fix = nee_get;
			return ret;
		}
		else {
			mainout << "Fail to change bounds, need new bounds to re-diagonalize =" << Ef << "\t" << "New Energy bound guesses: " << eelow << " -- " << eehigh << endl;
			nee_expect = min(int((eehigh - eelow) * averDos * 1.5), MaxSubspace - 10);
			ret = 1;
			return ret;
		}
	}else if(fmode == 1) {//use up to the first orbital that is occupied less than C_min (see syspara.hpp for def); here operations can be crude, real fix later;
		//find occupation Number and fermi energy
		double deltaNdiff; // this is the numerical difference between total occupied and the integer nee_want, should be small
		int icount = -1;//is the number of search in the fermi - distribution function to find the fermi energy		
		lastemin = ees[0]; nee_wantdecimal = nee_want;//here the particle number is fixed, so these values are known and given like this, these data not used, just for storage;
		//int FERMI_ENERGY_STATUS_NUM = findOccuCminLocateOneaboveEf(nee_want, nee_get, ees, eehigh, eelow, kBt, Ef, occu, deltaNdiff, icount, nee_fix, PotentialMin);
		
		int nee_expect_temp;
		int FERMI_ENERGY_STATUS_NUM = findOccuCminLocateOneaboveEf_tellExpect(nee_expect_temp, nee_want, nee_get, ees, eehigh, eelow, kBt, Ef, occu, deltaNdiff, icount, nee_fix, PotentialMin);
		if (FERMI_ENERGY_STATUS_NUM == 3 || FERMI_ENERGY_STATUS_NUM == 4) {
			nee_expect = min(int((eehigh - eelow) * averDos * 1.5), MaxSubspace - 10);//expect new Norb using existing dos
			nee_expect = max(nee_expect_temp, nee_expect);//yayun, 20240127, it happens that Nsubspace is too small when FERMI_ENERGY_STATUS_NUM == 4 "E:\sharedata01\yhu\workspace\data\CFDFTsquare\LXfracFixNee240116\28\slurm-35187_1.out"
			mainout << "Found: " << FERMI_ENERGY_STATUS_NUM << " fermiEnergy= " << Ef << "\t" << "New bound guesses: " << std::setprecision(6) << eelow << " -- " << eehigh << endl;									
			return ret;//next iteration
		}
		else if (FERMI_ENERGY_STATUS_NUM == 5) {
			nee_expect = min(int((eehigh - eelow) * averDos * 1.5), MaxSubspace - 10);//expect new Norb using existing dos
			mainout << "Special, large fermi when nee_want=nee_get, Found: FERMI_ENERGY_STATUS_NUM = " << FERMI_ENERGY_STATUS_NUM << " fermiEnergy=" << Ef << "\t" << "New bound guesses: " << std::setprecision(6) << eelow << " -- " << eehigh << endl;
			return ret;//next iteration
		}
		else if(FERMI_ENERGY_STATUS_NUM == 1){
			mainout << "Fail to find Ef, FERMI_ENERGY_STATUS_NUM = "<< FERMI_ENERGY_STATUS_NUM <<", fermiEnergy= " << Ef << ", icount is " << icount << ", deltaNdiff is " << std::setprecision(18) << deltaNdiff << "\t" << "New bound guesses: " << std::setprecision(6) << eelow << " -- " << eehigh << endl;
			nee_expect = min(int((eehigh - eelow) * averDos * 1.5), MaxSubspace - 10);//expect new Norb using existing dos	
			//nrun,nee_get,nee_expect,nee_fix,icount,deltaNdiff,eehigh,eelow,Ef
			fermiNotFound(filefermiNotFound,nrun,nee_get,nee_expect,nee_fix,icount,deltaNdiff,eehigh,eelow,Ef,FERMI_ENERGY_STATUS_NUM);
			ret = 1;
			return ret;//rediagonalize for fermi energy
		}
		else if (FERMI_ENERGY_STATUS_NUM == 2) {
			nee_expect = min(int((eehigh - eelow) * averDos * 1.5), MaxSubspace - 10);//expect new Norb using existing dos
			mainout << "Cmin not satisfied, FERMI_ENERGY_STATUS_NUM= " << FERMI_ENERGY_STATUS_NUM << " eigen low&high: " << ees[0] <<"--" << ees[nee_want - 1] <<  ",fermiEnergy= " << Ef << " New bound guesses: " << std::setprecision(6) << eelow << " -- " << eehigh << endl;
			fermiNotFound(filefermiNotFound, nrun, nee_get, nee_expect, nee_fix, icount, deltaNdiff, eehigh, eelow, Ef, FERMI_ENERGY_STATUS_NUM);
			ret = 1;
			return ret;//rediagonalize for a few more orbitals
		}		
	}
	else if (fmode == 2) {//find occupation Number
		int FERMI_ENERGY_STATUS_NUM = findOccuCminfixedEf(nee_want, nee_get, ees, eehigh, eelow, kBt, Ef, occu, nee_fix, PotentialMin, lastemin, nee_wantdecimal);
		if (FERMI_ENERGY_STATUS_NUM == 2) {//this condition is the only possibility
			nee_expect = min(int((eehigh - eelow) * averDos * 1.5), MaxSubspace - 10);//expect new Norb using existing dos
			mainout << std::setprecision(6) << " Emin: " << lastemin << "\t" << "New bound guesses: " << eelow << " -- " << eehigh << endl;
			return ret;//next iteration
		}
	}
	return ret;//will not happen
}

int calibrateFixEigen(int nee_want, int& nee_get, const double* ees, double& eehigh, double& eelow, double& kBt, double& Ef, double* occu, int& nee_expect, double& averDos, int& nfix) {
	using namespace std;
	using namespace MAIN;
	int ret = 0;
	double last;
	averDos = double(nee_get + 5.0) / (ees[nee_get - 1] - ees[0]); // maximum estimate of density of states using output energy

	//find occupation Number and fermi energy
	double deltaNdiff; int icount = -1;// this is the numerical difference between total occupied and the integer nee_want, should be small
	if (findOccupationNumber(nee_want, nee_get, ees, eehigh, eelow, kBt, Ef, occu, deltaNdiff, icount)) {
		//reportOccupation(-1, nee_get);
		mainout << "This is the temperature test, fermiEnergy=" << Ef << ",deltaNdiff is" << std::setprecision(18) << deltaNdiff << "\t" << "Old Energy bound guesses: " << std::setprecision(6) << eelow << " -- " << eehigh << endl;
	}
	else {
		mainout << "Fail to find fermi energy, fermiEnergy=" << Ef << ",deltaNdiff is" << std::setprecision(18) << deltaNdiff << "\t" << "Old Energy bound guesses: " << std::setprecision(6) << eelow << " -- " << eehigh << endl;
		eelow = ees[0] - CyclotronGap;
		eehigh = ees[nee_get - 1] + CyclotronGap;
		nee_expect = min(int((eehigh - eelow) * averDos * 1.5), MaxSubspace - 10);//expect new Norb using existing dos
		ret = 1;
		return ret;
	}
	//change the bounds. Also possibly change occupationNumber and fermiEnergy!!!!
	bool modifyFermiOccu = true;// This is to say when there are too many orbitals, then use only certain number of them to calculate ef and occu
	bool modified = false;//tell mainout if ef and occu modified
	if (shiftEboudFirstOccuTemperature(nee_want, nee_get, ees, eehigh, eelow, kBt, Ef, occu, modifyFermiOccu, modified)) {//shiftEboudFirstOccuTemperature
		mainout << "New energy bounds set=" << "\t" << "New Energy bound guesses: " << eelow << " -- " << eehigh << endl;
		if (modified) {
			//reportOccupation(-2, nee_get);
			mainout << "New occupations are reset as well as fermiEnergy, which is=" << Ef << endl;
		}
		nee_expect = min(int((eehigh - eelow) * averDos * 1.5), MaxSubspace - 10);
		return ret;
	}
	else {
		mainout << "Fail to change bounds, need new bounds to re-diagonalize =" << Ef << "\t" << "New Energy bound guesses: " << eelow << " -- " << eehigh << endl;
		nee_expect = min(int((eehigh - eelow) * averDos * 1.5), MaxSubspace - 10);
		ret = 1;
		return ret;
	}
}





/* Plain -laplacian. */
template<typename T>
void mLaplacian(T& ham){
	const auto ax2 = ax*ax, ay2=ay*ay;
	for(iint i=0; i<Nxy; i++) {
		const auto l=lnn(i), r=rnn(i), u=unn(i), d=dnn(i);
		if(l >= 0) {
			ham[At(l,i)] -= 1/ax2;
		}
		if(r >= 0) {
			ham[At(r,i)] -= 1/ax2;
		}
		if(u >= 0) {
			ham[At(u,i)] -= 1/ay2;
		}
		if(d >= 0) {
			ham[At(d,i)] -= 1/ay2;
		}
		ham[At(i,i)] -= (-2/ax2 -2/ay2);
	}
}


/* Covariant kinetic energy, aka (p-qA)^2/2m */
template<typename T>
void addKinetic(T& ham, const arma::mat &A_x, const arma::mat &A_y){
	const auto ax2 = ax*ax, ay2=ay*ay;
	iint ix, iy, jx, jy;
	for(iint i=0; i<Nxy; i++) {

		gridFromIndex(i, ix, iy);

		const auto l=lnn(i), r=rnn(i), u=unn(i), d=dnn(i); // indices of neighbors.

		/* multiply everything by the kinetic energy factor, since NN & host do not overlap */
		if(l >= 0) {
			gridFromIndex(l,jx,jy);
			ham[At(l,i)] += KE_scale*((-1/ax2)+I*(A_x(ix,iy)+A_x(jx,jy))/2.0/ax);//1/(2m^*) *hbar*eB_u/c [-d^2/dx^2 one of 4 components "-" -2IA_x d/dx one of 2 components]
		}
		if(r >= 0) {
			gridFromIndex(r,jx,jy);
			ham[At(r,i)] += KE_scale*((-1/ax2)-I*(A_x(ix,iy)+A_x(jx,jy))/2.0/ax);//1/(2m^*) *hbar*eB_u/c [-d^2/dx^2 one of 4 components "+" -2IA_x d/dx one of 2 components]
		}
		if(d >= 0) {
			gridFromIndex(d,jx,jy);
			ham[At(d,i)] += Mx_vs_My*KE_scale*((-1/ay2)+I*(A_y(ix,iy)+A_y(jx,jy))/2.0/ay);
		}
		if(u >= 0) {
			gridFromIndex(u,jx,jy);
			ham[At(u,i)] += Mx_vs_My*KE_scale*((-1/ay2)-I*(A_y(ix,iy)+A_y(jx,jy))/2.0/ay);
		}
		ham[At(i,i)] += KE_scale*((2/ax2+2/ay2) + pow(A_x(ix,iy),2) + pow(A_y(ix,iy),2));
	}
}

/* Covariant kinetic energy, aka (p-qA)^2/2m */
template<typename T>
void addKinetic_approx(T& ham, const arma::mat &A_x, const arma::mat &A_y){
	const auto ax2 = ax*ax, ay2=ay*ay;
	iint ix, iy, jx, jy;
	for(iint i=0; i<Nxy; i++) {

		gridFromIndex(i, ix, iy);

		const auto l=lnn(i), r=rnn(i), u=unn(i), d=dnn(i); // indices of neighbors.

		/* multiply everything by the kinetic energy factor, since NN & host do not overlap */
		if(l >= 0) {
			gridFromIndex(l,jx,jy);
			ham[At(l,i)] += KE_scale*((-1/ax2)+I*(A_x(ix,iy)+A_x(jx,jy))/2.0/ax);//1/(2m^*) *hbar*eB_u/c [-d^2/dx^2 one of 4 components "-" -2IA_x d/dx one of 2 components]
		}
		if(r >= 0) {
			gridFromIndex(r,jx,jy);
			ham[At(r,i)] += KE_scale*((-1/ax2)-I*(A_x(ix,iy)+A_x(jx,jy))/2.0/ax);//1/(2m^*) *hbar*eB_u/c [-d^2/dx^2 one of 4 components "+" -2IA_x d/dx one of 2 components]
		}
		if(d >= 0) {
			gridFromIndex(d,jx,jy);
			ham[At(d,i)] += Mx_vs_My*KE_scale*((-1/ay2)+I*(A_y(ix,iy)+A_y(jx,jy))/2.0/ay);
		}
		if(u >= 0) {
			gridFromIndex(u,jx,jy);
			ham[At(u,i)] += Mx_vs_My*KE_scale*((-1/ay2)-I*(A_y(ix,iy)+A_y(jx,jy))/2.0/ay);
		}
		ham[At(i,i)] += KE_scale*((2/ax2+2/ay2) + pow(A_x(ix,iy),2) + pow(A_y(ix,iy),2));
	}
}


/* Covariant kinetic energy, aka (p-qA)^2/2m */
template<typename T>
void addKinetic_Piers(T& ham, const arma::mat &A_x, const arma::mat &A_y){
	const auto ax2 = ax*ax, ay2=ay*ay;
	iint ix, iy, jx, jy;
	for(iint i=0; i<Nxy; i++) {
		/* (ix, iy) are row, c_(ix, iy)^dagger,  (jx, jy) are column, c_(jx, jy), ham elemet is t_[(ix, iy), (jx, jy)] c_(ix, iy)^dagger*c_(jx, jy)
		*/
		gridFromIndex(i, ix, iy);

		const auto l=lnn(i), r=rnn(i), u=unn(i), d=dnn(i); // indices of neighbors.

		/* multiply everything by the kinetic energy factor, since NN & host do not overlap */
		if(l >= 0) {
			gridFromIndex(l,jx,jy);
			ham[At(l,i)] += KE_scale*(-1/ax2)*exp(I*(A_x(ix,iy)+A_x(jx,jy))/2.0*ax);  
			/* because left, so the shift ax is negative, and there is also a negative sign in front of I
			*/
			//1/(2m^*) *hbar*eB_u/c [-d^2/dx^2 one of 4 components "-" -2IA_x d/dx one of 2 components]. IS a hopping when A==0.
		}
		if(r >= 0) {
			gridFromIndex(r,jx,jy);
			ham[At(r,i)] += KE_scale*(-1/ax2)*exp(-I*(A_x(ix,iy)+A_x(jx,jy))/2.0*ax);//1/(2m^*) *hbar*eB_u/c [-d^2/dx^2 one of 4 components "+" -2IA_x d/dx one of 2 components]
		}
		if(d >= 0) {
			gridFromIndex(d,jx,jy);
			ham[At(d,i)] += Mx_vs_My*KE_scale*(-1/ay2)*exp(I*(A_y(ix,iy)+A_y(jx,jy))/2.0*ay);
		}
		if(u >= 0) {
			gridFromIndex(u,jx,jy);
			ham[At(u,i)] += Mx_vs_My*KE_scale*(-1/ay2)*exp(-I*(A_y(ix,iy)+A_y(jx,jy))/2.0*ay);
		}
		//ham[At(i,i)] += KE_scale*((2/ax2+2/ay2) + pow(A_x(ix,iy),2) + pow(A_y(ix,iy),2));
		ham[At(i,i)] += KE_scale*(2/ax2)+Mx_vs_My*KE_scale*(2/ay2);
	}
}

/* Covariant kinetic energy, aka (p-qA)^2/2m */
template<typename T>
void addKinetic_Piers_negative(T& ham, const arma::mat &A_x, const arma::mat &A_y){
	const auto ax2 = ax*ax, ay2=ay*ay;
	iint ix, iy, jx, jy;
	for(iint i=0; i<Nxy; i++) {
		/* (ix, iy) are row, c_(ix, iy)^dagger,  (jx, jy) are column, c_(jx, jy), ham elemet is t_[(ix, iy), (jx, jy)] c_(ix, iy)^dagger*c_(jx, jy)
		*/
		gridFromIndex(i, ix, iy);

		const auto l=lnn(i), r=rnn(i), u=unn(i), d=dnn(i); // indices of neighbors.

		/* multiply everything by the kinetic energy factor, since NN & host do not overlap */
		if(l >= 0) {
			gridFromIndex(l,jx,jy);
			ham[At(l,i)] += KE_scale*(-1/ax2)*exp(-I*(A_x(ix,iy)+A_x(jx,jy))/2.0*ax);  
			/* because left, so the shift ax is negative, and there is also a negative sign in front of I
			*/
			//1/(2m^*) *hbar*eB_u/c [-d^2/dx^2 one of 4 components "-" -2IA_x d/dx one of 2 components]. IS a hopping when A==0.
		}
		if(r >= 0) {
			gridFromIndex(r,jx,jy);
			ham[At(r,i)] += KE_scale*(-1/ax2)*exp(I*(A_x(ix,iy)+A_x(jx,jy))/2.0*ax);//1/(2m^*) *hbar*eB_u/c [-d^2/dx^2 one of 4 components "+" -2IA_x d/dx one of 2 components]
		}
		if(d >= 0) {
			gridFromIndex(d,jx,jy);
			ham[At(d,i)] += Mx_vs_My*KE_scale*(-1/ay2)*exp(-I*(A_y(ix,iy)+A_y(jx,jy))/2.0*ay);
		}
		if(u >= 0) {
			gridFromIndex(u,jx,jy);
			ham[At(u,i)] += Mx_vs_My*KE_scale*(-1/ay2)*exp(I*(A_y(ix,iy)+A_y(jx,jy))/2.0*ay);
		}
		//ham[At(i,i)] += KE_scale*((2/ax2+2/ay2) + pow(A_x(ix,iy),2) + pow(A_y(ix,iy),2));
		ham[At(i,i)] += KE_scale*(2/ax2)+Mx_vs_My*KE_scale*(2/ay2);
	}
}

/* Covariant kinetic energy, aka (p-qA)^2/2m */
template<typename T>
void addPeriodKineticAndOnsite(T& ham, const arma::mat &A_x, const arma::mat &A_y, const bool usePertemp){
	if (usePertemp == true){
		const auto ax2 = ax*ax, ay2=ay*ay;
		iint ix, iy, jx, jy;
		for(iint i=0; i<Nxy; i++) {

			gridFromIndex(i, ix, iy);

			if (ix == 0){
				const auto l = gridToIndex(Nx -1, iy); //left neighbor index
				gridFromIndex(l,jx,jy);
				ham[At(l,i)] += KE_scale*((-1/ax2)+I*(A_x(ix,iy)+A_x(jx,jy))/2.0/ax);//1/(2m^*) *hbar*eB_u/c [-d^2/dx^2 one of 4 components "-" -2IA_x d/dx one of 2 components]			
			}
			else if(ix == Nx-1){
				const auto r = gridToIndex(0, iy); //right neighbor index
				gridFromIndex(r,jx,jy);
				ham[At(r,i)] += KE_scale*((-1/ax2)-I*(A_x(ix,iy)+A_x(jx,jy))/2.0/ax);//1/(2m^*) *hbar*eB_u/c [-d^2/dx^2 one of 4 components "+" -2IA_x d/dx one of 2 components]

			}

		}
	}
}




/* Add generic onsite potential V to hamiltonian matrix ham */
template<typename T>
void addOnsite(T& ham, const arma::mat &V){
	iint ix, iy;
	for(iint i=0; i<Nxy; i++) {
		gridFromIndex(i, ix, iy);
		ham[At(i,i)] += V(ix,iy);
	}
}

/* Add Ion Confinement energy. */
template<typename T>
void addIonConfinement(T& ham, const arma::mat &V){
	iint ix, iy;
	for(iint i=0; i<Nxy; i++) {
		gridFromIndex(i, ix, iy);
		ham[At(i,i)] += V(ix,iy)*scaleIon;
	}
}


/* Add Hartree potential. */
template<typename T>
void addHartree(T& ham, const arma::mat &V){
	iint ix, iy;
	for(iint i=0; i<Nxy; i++) {
		gridFromIndex(i, ix, iy);
		ham[At(i,i)] += V(ix,iy)*scaleHartree;
	}
}
/*
 * Calculates total energy, using states in Eigenvecs in the meanfield density seaDensity.
 * Uses and resets the hamiltonian in hamop.
 * with temperature. This condiers the Exc with filling factor calculated in real external B
 */
double calcEnergyTemperatureRealB(double& Etotal, double& KE, double& Eext, double& EHartree, double& Excsum, const MMSpmatR_op<iint>& hamop, 
	const arma::mat& seaDensity, const arma::mat& ionCoulomb, const arma::cx_mat& Eigenvecs, const iint nee_get, const double* occu, const double& deltaB) {
	using namespace arma;

	// reset ham
	// std::fill(hamop.spmat->vals, hamop.spmat->vals + hamop.spmat->nnz, (ComplexF64)(0));//yayun,20240107,do not change hamop because will saveHamToFile
	/*I decide not to update these potentials for calculating energy, YY210416;
	calcA(seaDensity, deltaB); 
	AddCircumcircleA(useCircumcircleA, deltaB);
	calcCoulomb(seaDensity);
	*/
	//addKinetic(*(hamop.spmat), VECTORPOT::A_x_mat, VECTORPOT::A_y_mat); //yayun,20240107,do not change hamop because will saveHamToFile
	//addPeriodKineticAndOnsite(*(hamop.spmat), VECTORPOT::A_x_mat, VECTORPOT::A_y_mat, shape_fill_num);
	//auto H_ev = cx_mat(Nxy, nee_get);//yayun,20240107,do not change hamop because will saveHamToFile
	//mkl_mm<iint>(hamop, nee_get, reinterpret_cast<const MKL_Complex16*>(Eigenvecs.memptr()), reinterpret_cast<MKL_Complex16*>(H_ev.memptr()));// is this matrix more efficient than one by one kinetic energy?
	//
	KE = 0.0;
	/*
	for (iint i = 0; i < nee_get - 1; i++) {
		KE += real(cdot(Eigenvecs.col(i), H_ev.col(i))) * occu[i];
	}	
	*///yayun,20240107,do not change hamop because will saveHamToFile
	//KE = real(cdot(Eigenvecs.cols(0, nee_get - 1), H_ev));
	EHartree = accu(seaDensity % COULOMB::VCoulomb / 2.0) ; // % means .*, * scaleHartree not multiplied because VCoulomb now is used as it is;
	Eext = accu(seaDensity % ionCoulomb); // % means .*, not multiplying  * scaleIon because ionCoulomb is used as it is;
	Excsum = 0.0;
	seaDensity.for_each([&Excsum, &deltaB](const double& val) {Excsum += ExcRealFilling(val, deltaB); });
	Excsum *= scaleXC;
	Etotal = KE + EHartree + Eext + Excsum;

	return Etotal;
}

//210411 commented out, Yayun; 
arma::mat XX(Nx, Ny), YY(Nx, Ny); // for debugging  ?????what is this? This is global quantity


/* Calculates p-qA the dynamical momentum, finite temperature with occupationNumber */
void calcDynMomentumTemperature(arma::mat& dynP_x, arma::mat& dynP_y, const arma::mat& rho, const arma::cx_mat& evs, const int nee_get, const double* occu, const arma::mat& A_x, const arma::mat& A_y) {

	using namespace arma;

	for (iint i = 0; i < Nx; i++) {
		for (iint j = 0; j < Ny; j++) {
			iint ii = gridToIndex(i, j); // index along eigenvector, p-qA at each position ii

			/* Calculate A * psi^2 */
			const auto APsi_x = rho(i, j) * A_x(i, j);
			const auto APsi_y = rho(i, j) * A_y(i, j);


			/* Calculate \sum_{k=0}^{np-1} psi_k^*(r)|P|psi_k(r) */
			const auto l = lnn(ii), r = rnn(ii), u = unn(ii), d = dnn(ii); // indices of neighbors
			ComplexF64 PxPsi = 0, PyPsi = 0;
			for (iint k = 0; k < nee_get; k++) {// sum over all occupied orbitals
				if (l >= 0) {
					PxPsi += conj(evs(ii, k)) * -I * (evs(ii, k) - evs(l, k)) / 2.0 / ax * occu[k];
				}
				if (r >= 0) {
					PxPsi += conj(evs(ii, k)) * -I * (evs(r, k) - evs(ii, k)) / 2.0 / ax * occu[k];
				}
				if (d >= 0) {
					PyPsi += conj(evs(ii, k)) * -I * (evs(ii, k) - evs(d, k)) / 2.0 / ay * occu[k];
				}
				if (u >= 0) {
					PyPsi += conj(evs(ii, k)) * -I * (evs(u, k) - evs(ii, k)) / 2.0 / ay * occu[k];
				}
			}

			XX(i, j) = imag(PxPsi);
			YY(i, j) = imag(PyPsi);

			dynP_x(i, j) = APsi_x + real(PxPsi);
			dynP_y(i, j) = APsi_y + real(PyPsi);
		}
	}
}



/* Calculates p-qA the dynamical momentum. */
void calcDynMomentum(arma::mat &dynP_x, arma::mat &dynP_y, const arma::mat &rho, const arma::cx_mat &evs, const iint np, const arma::mat &A_x, const arma::mat &A_y) {

	using namespace arma;

	for(iint i = 0; i < Nx; i++) {
		for(iint j = 0; j < Ny; j++) {
			iint ii = gridToIndex(i,j); // index along eigenvector

			/* Calculate A * psi^2 */
			const auto APsi_x = rho(i,j)*A_x(i,j);
			const auto APsi_y = rho(i,j)*A_y(i,j);


			/* Calculate \sum_{k=0}^{np-1} psi_k^*(r)|P|psi_k(r) */
			const auto l=lnn(ii), r=rnn(ii), u=unn(ii), d=dnn(ii); // indices of neighbors
			ComplexF64 PxPsi = 0, PyPsi = 0;
			for(iint k = 0; k < np; k++) {
				if(l >= 0) {
					PxPsi += conj(evs(ii,k))*-I*(evs(ii,k)-evs(l,k))/2.0/ax;
				}
				if(r >= 0) {
					PxPsi += conj(evs(ii,k))*-I*(evs(r,k)-evs(ii,k))/2.0/ax;
				}
				if(d >= 0) {
					PyPsi += conj(evs(ii,k))*-I*(evs(ii,k)-evs(d,k))/2.0/ay;
				}
				if(u >= 0) {
					PyPsi += conj(evs(ii,k))*-I*(evs(u,k)-evs(ii,k))/2.0/ay;
				}
			}

			XX(i,j) = imag(PxPsi);
			YY(i,j) = imag(PyPsi);

			dynP_x(i,j) = APsi_x+real(PxPsi);
			dynP_y(i,j) = APsi_y+real(PyPsi);
		}
	}
}


#endif
