// This file is under public domain.
#ifndef MAIN_UTIL_HPP
#define MAIN_UTIL_HPP

#include<algorithm> // for copy

template<typename T>
void resizePointerArray(T* a, const size_t n, const size_t  new_n)//copy the first n elements of a into a resized a of size new_n;n<=new_n
{
	T* new_a = new T[new_n];
	copy(a, a + n, new_a);
	delete[] a;
	a = new_a;
}

void writePotentials() {
	using namespace arma;
	using namespace VECTORPOT;
	{
		ofstream f(prefix + "greeAx");

		for (auto j = 0; j < Ny; j++) {

			for (auto i = 0; i < Nx; i++) {
				f << VECTORPOT::greenA_x(i, j) << '\t';
			}
			f << endl;
		}
	}
	{
		ofstream f(prefix + "A_x_mat");

		for (auto j = 0; j < Ny; j++) {
			for (auto i = 0; i < Nx; i++) {
				f << A_x_mat(i, j) << '\t';
			}
			f << endl;
		}
	}
	{
		ofstream f(prefix + "A_y_mat");

		for (auto j = 0; j < Ny; j++) {
			for (auto i = 0; i < Nx; i++) {
				f << A_y_mat(i, j) << '\t';
			}
			f << endl;
		}
	}
	{
		ofstream f(prefix + "Vcoulomb");
		
		for (auto j = 0; j < Ny; j++) {
			for (auto i = 0; i < Nx; i++) {
				f << COULOMB::VCoulomb(i, j) << '\t';
			}
			f << endl;
		}
	}
	
}

void writeHam() {
	using namespace arma;
	INFO(Nx);
	INFO(Ny);
	INFO(Nxy);
	INFO(ham.nnz);
	{
		ofstream f(prefix + "mat");
		f << std::setprecision(18);
		for (auto i = 0; i < ham.nnz; i++) {
			f << ham.rows1[i] << '\t' << ham.cols1[i] << '\t' << real(ham.vals[i]) << '\t'
			  << imag(ham.vals[i]) << endl;
		}
	}
	{
		ofstream f(prefix + "mat.ana");
		f << std::setprecision(18);
		f << "#row\tcol\tv.real\tv.imag\tix\tiy\tjx\tjy" << endl;// This is the header line or format of writtings
		for (auto i = 0; i < ham.nnz; i++) {
			int ix, iy, jx, jy;
			gridFromIndex(ham.rows[i], ix, iy);
			gridFromIndex(ham.cols[i], jx, jy);
			f << ham.rows1[i] << '\t' << ham.cols1[i] << '\t' << real(ham.vals[i]) << '\t'
			  << imag(ham.vals[i]) << "\t" << ix << '\t' << iy << '\t' << jx << '\t'<<jy<<endl;
		}
	}
}

void printHaminfo() {
	using namespace arma;
	for(int i=0; i<Nxy+2; i++) {
		MAIN::mainout << ham.rowstart[i] << endl;
	}
	MAIN::mainout << ham.n_rows << ' ' << ham.n_cols << endl;
	for(auto i = 0; i < ham.nnz; i++) {
		MAIN::mainout << ham.rows[i] <<' '<< ham.cols[i] <<' '<< ham.vals[i] << '\n';
		MAIN::mainout << endl;
	}
}

void reporteigen(int n, int nee)// both eigenvalues and eigenvectors
{
	using namespace std;
	using namespace arma;
	string const sseq = (n == -1) ? "lastest" : std::to_string((long long)n);
	{
		// write eigenval
		ofstream f(prefix + "ees_" + sseq + ".txt");

		for (auto i = 0; i < nee; i++)
		{
			f << Eigenvals[i] << endl;
		}
	}

	if(n>0){
		// write eigenvecs
		cx_mat a = Eigenvecs.cols(0, nee - 1);

		string fname = prefix + "evs_" + sseq + ".bin";

		writenumarray(fname, a);
	}
}

void reporteigenvec(int n, int nee)
{
	using namespace std;
	using namespace arma;
	string const sseq = std::to_string((long long)n);
		// write eigenvecs
		cx_mat a = Eigenvecs.cols(0, nee - 1);

		string fname = prefix + "evs_" + sseq + ".bin";

		writenumarray(fname, a);
}

void reporteigenvecFinished(const string& sseq,const int& nee)
{
	using namespace std;
	using namespace arma;
	// write eigenvecs
	cx_mat a = Eigenvecs.cols(0, nee - 1);
	writenumarray(sseq, a);
}

void reporteigenvalues(int n, int nee)
{
	using namespace std;
	using namespace arma;
	string const sseq = std::to_string((long long)n);
	{
		// write eigenval
		ofstream f(prefix + "ees_" + sseq + ".txt");

		for (auto i = 0; i < nee; i++)
		{
			f << Eigenvals[i] << endl;
		}
	}
}

void reporteigenvaluesFeast(int n, int nee)//
{
	using namespace std;
	using namespace arma;
	string const sseq = std::to_string((long long)n);
	{
		// write eigenval
		ofstream f(prefix + "Feast_" + sseq + ".txt");

		for (auto i = 0; i < nee; i++)
		{
			f << Eigenvals[i] << endl;
		}
	}
}

void reportOccupation(int n, int nee) { 
	using namespace std;
	using namespace arma;
	string const sseq = std::to_string((long long)n);
	{
		// write eigenval
		ofstream f(prefix + "occupation" + sseq + ".txt");

		for (auto i = 0;i < nee;i++) {
			f << occupationNumber[i] << endl;
		}
	}

}

void occupationDiff(const int &n, const int &nee, double &x, arma::rowvec& y) {
	using namespace std;
	for (auto i = 0; i < n; i++) {
		 y.col(i)=occupationNumber[i] - occupationNumberOld[i];
		//update eig and occu Old files
	}
	x = accu(abs(y.cols(0,n-1)));
	std::copy(occupationNumber, occupationNumber + nee, occupationNumberOld);
}


void eesDiff(const int &n, const int &nee, double& x, arma::rowvec& y) {
	using namespace std;
	for (int i = 0; i < n; i++) {
		y.col(i) = Eigenvals[i] - EigenvalsOld[i];
		//update eig and occu Old files
	}
	x = accu(abs(y.cols(0, n - 1)));
	std::copy(Eigenvals, Eigenvals + nee, EigenvalsOld);
}

void evsOverlapDiff(const int& n, const int& nee, arma::rowvec & x, arma::cx_rowvec& y, double& z) {/*absolute diff, overlap*/
	using namespace std;
	using namespace arma;
	arma::cx_mat phasemat(n, n);
	arma::cx_rowvec phase(1, n);
	urowvec phaseindex(1, n);
	phasemat = trans(EigenvecsOld.cols(0, n - 1)) * Eigenvecs.cols(0, n - 1);//expansion coefficient in old eigenstate
	phase = max(phasemat,0);
	y.cols(0, n - 1) = phase;
	phaseindex = index_max(phasemat, 0);//pick out the largest-weighted old column
	phase = phase/ abs(phase); //the coefficient should be considered either only phase or in full amplitude (not divided by abs)
/*get rid of phase diff*/
	EigenvecsOld.cols(0, n - 1) = EigenvecsOld.cols(phaseindex);//ready to compare with the largest-weighted old column
	EigenvecsOld.cols(0, n - 1).each_row() %= phase;//try to take care of trivial phase difference
	EigenvecsOld.cols(0, n - 1) = Eigenvecs.cols(0, n - 1) - EigenvecsOld.cols(0, n - 1);// bare difference, not considering phase here if above line is removed
	x.cols(0, n - 1)=sum(pow(abs(EigenvecsOld.cols(0, n - 1)), 2), 0);
	z = accu(x.cols(0, n - 1));
	EigenvecsOld.resize(Nxy, nee);//yayun 240219, to save storage such to avoid segmentation fault
	EigenvecsOld.cols(0, nee - 1) = Eigenvecs.cols(0,nee-1);
}

void newUpdateRate(bool& ztoolargelast, int& npre, const int& n, double& up, double& rate, double& theta, double& expansion, double& delta,
	arma::rowvec& x, arma::rowvec& y/*occudelta*/, const arma::rowvec& z/*evsdiff*/, const arma::rowvec& u/*occupationNumber*/, const double& VTdiff, 
	arma::rowvec& v, arma::rowvec& w, double& rateE, double& thetaE, double& expansionE, double& deltaE) {
	using namespace std;
	bool ztoolarge = max(z.cols(0, n - 1) % u.cols(0, n - 1)) > 0.05;//This might be still too conservative
	//MAIN::mainout << "ztoolarge" << max(z.cols(0, n - 1) % u.cols(0, n - 1)) << endl;
	//MAIN::mainout << "ztoolarge" << ztoolarge << endl;
	double Theta1 = Pi / 2*0.25;
	npre = min(npre, n);

	//ees angle
	arma::rowvec vweight=v.cols(0, npre - 1) % u.cols(0, npre - 1);
	arma::rowvec wweight = w.cols(0, npre - 1) % u.cols(0, npre - 1);	
	double lenv = sqrt(accu(pow(vweight.cols(0, npre - 1), 2)));//pre
	double lenw = sqrt(accu(pow(wweight.cols(0, npre - 1), 2)));//this
	theta = dot(vweight.cols(0, npre - 1), wweight.cols(0, npre - 1));
	expansionE = theta / pow(lenv, 2);
	theta = theta / lenv / lenw;
	if (theta<1.0 && theta>-1.0) {
		theta = acos(theta);
	}
	else if (theta >= 1.0 && theta <= 1.0001) {
		theta = 0.0;
	}
	else if (theta <= -1.0 && theta > -1.0001) {
		theta = Pi;
	}
	else {
		MAIN::mainout << "wrong value of thetaE" << theta << endl;
	} 
	delta = sqrt(accu(pow(vweight.cols(0, npre - 1) - wweight.cols(0, npre - 1), 2))) / lenv;// difference between the two changes relative to the previous change, all measured by magnitude
	if (theta < Theta1 && !ztoolarge) {
		rate = 0.3 * pow(sin((1 + theta / Theta1) * Pi / 2), 6) * (tanh(-delta) + 1);//slow increase, pow 6
	}
	else if (theta >= Theta1) {
		if (up > UpdateRatioBase) {
			//up = UpdateRatioBase;
			up *= 0.5;
		}
		rate = -0.2 * pow(sin((1 + (Pi - theta) / (Pi - Theta1)) * Pi / 2), 4) * tanh(lenw / lenv);//faster decrease, pow 4
	}
	up *= (1 + rate);
	//hard UpdateRatioLowBound, consider ees sudden change??
	if (up > UpdateRatioUpperBound) {
		if (delta > 1.0) {
			up = max(UpdateRatioLowBound,up / delta);//UpdateRatioUpperBound;
		}
	}
	rateE = rate;
	thetaE = theta;
	deltaE = delta;

	Theta1 = Pi / 2 * 0.25;//Theta1 for occu angle
	double lenx = sqrt(accu(pow(x.cols(0, npre - 1), 2)));//pre
	double leny = sqrt(accu(pow(y.cols(0, npre - 1), 2)));//this
	theta = dot(x.cols(0, npre - 1), y.cols(0, npre - 1));
	expansion = theta / pow(lenx, 2);
	theta = theta / lenx / leny;
	if (theta<1.0 && theta>-1.0) {
		theta = acos(theta);
	}
	else if (theta >= 1.0 && theta <= 1.0001) {
		theta = 0.0;
	}
	else if (theta <= -1.0 && theta > -1.0001) {
		theta = Pi;
	}
	else {
		MAIN::mainout << "wrong value of theta" << theta << endl;
	} 
	delta = sqrt(accu(pow(x.cols(0, npre - 1) - y.cols(0, npre - 1),2))) / lenx;// difference between the two changes relative to the previous change, all measured by magnitude
	if (theta < Theta1 && !ztoolarge) {
		rate=0.3 * pow(sin((1 + theta / Theta1) * Pi / 2), 4) * (tanh(-delta) + 1);//slow increase, pow 6 now 4
	}else if(theta >= Theta1) {
		if (up> UpdateRatioBase){
		//up = UpdateRatioBase;
		up *= 0.5;
		}
	    rate= -0.2 * pow(sin((1 + (Pi - theta) / (Pi - Theta1)) * Pi / 2),4) * tanh(leny/lenx);//faster decrease, pow 4
	}
	up *= (1 + rate);
	if ((ztoolarge  || ztoolargelast) && up > UpdateRatioBase) {//; //This is when single orbitals are not stable
		if (theta >= Theta1) {
			up *= 0.4;
		}
		else {
			up *= 0.7;// many changes in single particle orbitals occur, maybe just trivially because system rotating and orbitals too degenerate
		}
		//MAIN::mainout << "ztoolarge up in" << endl;
	}

	//now consider absolute change of occupation for the following two applications
	lenx = accu(abs(x.cols(0, npre - 1)));//pre, reuse of varibales
	leny = accu(abs(y.cols(0, n - 1)));//this
	lenx = max(lenx, leny);//take the maximum one
	//also consider absolute change of wave fucntions for the following two applications
	leny = accu(z.cols(0, n - 1) % u.cols(0, n - 1));//wave function changes weighted sum, not the same as but similar to vectorDiff
	lenx = max(lenx, leny);//take the maximum one again, though not the equivalent thing; did not further consider ees change here because do not know how to decide if it is large

	//if (lenx>0.5) {// change in occupation larger than 0.01 is too much
	//	up *= 0.5;
	//}
	if (up >0.8 ) {//&& VTdiff>0.01, no, do not let mixing exceed 1 //This is to consider mixing larger than 1, which may in particular add into VT too much unnecessary new components. 
		//Accept mixing>1 only when VT do not change much relatively
		up = 0.8;// what if too fast???!
	}

	if (up < UpdateRatioLowBound ) {//This is to spefically increase mixing when absolute change in occu is small while mixng is minimized by theta=pi
		//up = UpdateRatioBase;
		up = UpdateRatioLowBound;// what if mixing too slow, may happen at early stage of iteration when max(lenx,leny) is large, but will deal with the late stage in the following when max(lenx,leny) is small
		if (lenx < pow(0.1,3) && VTdiff< pow(0.1, 4)) {// This is to avoid waste of time when occu is not going anywhere
			up = UpdateRatioBase;// change in occupation is so small, that I prefer to move on rapidly and see where things go; This is because here mainly occurs when there is always a theta=pi 
			if ( lenx < pow(0.1, 7) && VTdiff < pow(0.1, 5)) {//This can happen when use high accuracy in occuMinimumFirst
				up = 0.4;//Though theta is negative, the difference is small anyway, no big deal
				if ( lenx < pow(0.1, 11) && VTdiff < pow(0.1, 6)) {//This can happen when use high accuracy in occuMinimumFirst
					up = 0.8;//Though theta is negative, the difference is small anyway, no big deal
				}
			}
		}
	}

	//hard UpdateRatioLowBound
	if (up > UpdateRatioUpperBound) {
		if (delta>1.0) {
			up = max(UpdateRatioLowBound, up / delta);
			//up = UpdateRatioUpperBound;
		}
	}

	//up = 0.05;//do nothing

	x.cols(0, n - 1) = y.cols(0, n - 1);
	v.cols(0, n - 1) = w.cols(0, n - 1);
	npre = n;
	ztoolargelast = ztoolarge;
}


/* write XC potential. */
void writeVxc(int n, const arma::mat& val) {
	using namespace std;
	using namespace arma;
	arma::mat Vxctemp(size(val));
	for (auto j = 0; j < Ny; j++) {
		for (auto i = 0; i < Nx; i++) {
			Vxctemp(i, j) = Vxc(val(i, j)) * scaleXC;
		}
	}
	Vxctemp.save(prefix + "Vxc" + to_string((long long)n) + ".csv", arma::csv_ascii);
}

/* Record a row of energy entries of the current iteration. */
void recordEnergy(const string& filename, const int& nrun, const double& Etotal, const double& KE, const double& Eext, const double& EHartree, const double& Excsum) {
	using namespace std;
	ofstream f(filename, nrun == 0 ? ios_base::out : (ios_base::out | ios_base::app));
	if (nrun == 0) {
		f << "# nrun\tEtotal\tKE\tEext\tEHartree\tExcsum" << endl;
	}
	f << std::setprecision(18);
	f << nrun << '\t' << Etotal << '\t' << KE << '\t' << Eext
		<< '\t' << EHartree << '\t' << Excsum << endl;
	f.close();
}

/* Record a row of status including energy entries etc of the current iteration. */
void recordfermi(const string& filename, const int& nrun, const double& a, const double& b, const double& c, const double& d, const double& e,
	const double& f, const double& g, const double& h, const double& i, const double& j, const double& k, const double& l, const double& m,
	const double& n, const double& o, const double& p, const double& q, const double& r, const double& s, const double& t, const double& u, const double& v, const double& w,
	const double& x) {
	using namespace std;
	ofstream file(filename, nrun == 0 ? ios_base::out : (ios_base::out | ios_base::app));

	if (nrun == 0) {
		file << "# nrun\tNee\tNsubspace\tminNee\tEnergyhigh\tEnergylow\tfermiEnergy\toccudiff\teesdiff\trelVTDiff\trelVTDiffAim\trelDensityDiff" << 
			"\trelDensityDiffAim\tvectorDiff\tUpdateRatio\trateOccu\tthetaOccu\texpansionOccu\tdeltaOccu\tBigPeriodDiff\teesrate\teestheta\teesexpansion" << 
			"\teesdelta\ttotalPotentialMinimum" << endl;
	}
	file << std::setprecision(18);
	file << nrun << '\t' << a << '\t' << b << '\t' << c
		<< '\t' << d << '\t' << e << '\t' << f << '\t' << g << '\t' << h << '\t' << i << '\t' << j << '\t' << k << '\t' << l << '\t' << m <<
		'\t'<< n << '\t' << o << '\t' << p << '\t' << q << '\t' << r << '\t' << s << '\t' << t << '\t' << u << '\t' << v << '\t' << w << '\t' << x <<endl;
	file.close();
}

/* Record if fermi goes wrong iteration. */
void fermiNotFound(const string& filename, const int& nrun, const int& a, const int& b, const int& c, const int& d, const double& e, const double& f, const double& g, const double& h, const int& i) {
	using namespace std;
	ofstream file(filename, (ios_base::out | ios_base::app));
	//fermiNotFound(nrun, nee_get, nee_expect, nee_fix, icount, deltaNdiff, eehigh, eelow, Ef);
		file << std::setprecision(10);
		file << nrun << '\t' << a << '\t' << b << '\t' << c
			<< '\t' << d << '\t' << e << '\t' << f << '\t' << g << '\t' << h << '\t' << i << endl;
	file.close();
}

/* Record eigenvectors diff. */
void recordEigendiff(const string& filename, const int& nrun, const int& n, const double* x) {//const arma::rowvec&
	using namespace std;
	using namespace arma;
	ofstream file(filename, nrun == 1 ? ios_base::out : (ios_base::out | ios_base::app));
	file << std::setprecision(18);
	file << nrun << ',' << n;
		//file << '\t' << x.cols(0,n-1);
	for (int i = 0; i < n; i++) {
		file << ',' << x[i];
	}
	//mainout << x[1] << "double"<< endl;
	file << endl;
	file.close();
}


/* Record eigenvectors diff in binary. */
void recordEigendiffbinary(const string& filename, const int& nrun, const int& n, const arma::rowvec& x) {
	using namespace std;
	ofstream file(filename, nrun == 1 ? (ios_base::out | ios::binary) : (ios_base::out | ios_base::app | ios::binary));
	auto ptr = x.memptr();
	size_t len = sizeof(ptr[0]) * n;
	file.write((char*)ptr, len);
	file.close();
}

/* Record eigenvectors diff in binary. */
void recordEigenoverlapbinary(const string& filename, const int& nrun, const int& n, const arma::cx_rowvec& x) {
	using namespace std;
	ofstream file(filename, nrun == 1 ? (ios_base::out | ios::binary) : (ios_base::out | ios_base::app | ios::binary));
	auto ptr = x.memptr();
	size_t len = sizeof(ptr[0]) * n;
	file.write((char*)ptr, len);
	file.close();
}

/* Record eigenvectors overlap. */
void recordEigenoverlap(const string& filename, const int& nrun, const int& n, complex<double>* x) {//const arma::cx_rowvec& x
	using namespace std;
	using namespace arma;
	ofstream file(filename, nrun == 1 ? ios_base::out : (ios_base::out | ios_base::app));
	file << std::setprecision(18);
	file << nrun << ',' << n;
	//mainout << x[1] << "complex" <<endl;
	//mainout << x[2] << "complex" << endl;
	for (int i = 0; i < n; i++) {
		file << ',' << real(x[i])<<'+'<<imag(x[i]) <<'j';
	}
	file << endl;
	file.close();
}


/* Record occupations of converged/finished table parameters in one file. */
void recordoccuConvergeOrFinished(const string& filename, const int& tableNumtemp, const int& nee) {
	using namespace std;
	ofstream file(filename, tableNumtemp == 0 ? ios_base::out : (ios_base::out | ios_base::app));
	file << std::setprecision(18);
	file << tableNumtemp << ',' << nee;
	for (int i = 0; i < nee; i++) {
		file << ',' << occupationNumber[i];
	}
	file << endl;
	file.close();
}

/* Record ees of converged/finished table parameters in one file. */
void recordEigenConvergeOrFinished(const string& filename, const int& tableNumtemp, const int& nee) {
	using namespace std;
	ofstream file(filename, tableNumtemp == 0 ? ios_base::out : (ios_base::out | ios_base::app));
	file << std::setprecision(18);
	file << tableNumtemp << ',' << nee;
	for (int i = 0; i < nee; i++) {
		file << ',' << Eigenvals[i];
	}
	file << endl;
	file.close();
}



/* Record occupations in one file. */
void recordoccuAll(const string& filename, const int& nrun, const int& nee) {
	using namespace std;
	ofstream file(filename, nrun == 0 ? ios_base::out : (ios_base::out | ios_base::app ));
	file << std::setprecision(18);
	file << nrun << ',' << nee;
	for (int i = 0; i < nee; i++) {
		file << ',' << occupationNumber[i];
	}
	file << endl;
	file.close();
}

/* Record occupations in one file in binary. */
void recordoccuAllbinary(const string& filename, const int& nrun, const int& nee) {
	using namespace std;
	ofstream file(filename, nrun == 0 ? (ios_base::out | ios::binary) : (ios_base::out | ios_base::app | ios::binary));
	auto ptr = occupationNumber;
	size_t len = sizeof(ptr[0]) * nee;
	file.write((char*)ptr, len);
	file.close();
}

/* Record ees in one file. */
void recordEigenAll(const string& filename, const int& nrun, const int& nee) {
	using namespace std;
	ofstream file(filename, nrun == 0 ? ios_base::out : (ios_base::out | ios_base::app));
	file << std::setprecision(18);
	file << nrun << ',' << nee;
	for (int i = 0; i < nee; i++) {
		file << ',' << Eigenvals[i];
	}
	file << endl;
	file.close();
}

/* Record ees in one file in binary. */
void recordEigenAllbinary(const string& filename, const int& nrun, const int& nee) {
	using namespace std;
	ofstream file(filename, nrun == 0 ? (ios_base::out | ios::binary) : (ios_base::out | ios_base::app | ios::binary));
	auto ptr = Eigenvals;
	size_t len = sizeof(ptr[0]) * nee;
	file.write((char*)ptr, len);
	file.close();
}

//loading old results when starting a new tableNum that has not been started before, only want to make sure I start with previous density and VT
/*loadPreviousTableResults almost same as loadPreviousRoundResults, apart from the reference parameter name: tableNumnew, and its usage*/
bool loadPreviousTableResults(const int& tableNumnew, arma::mat& preTableSeaDensity, arma::mat& preTableVT, double& energylowtemp, double& energyhightemp, int& Nsubspacetemp,
	const arma::umat& TableOrdertempInversetemp, const arma::umat& TableOrderUseOldResultstemp, const arma::mat& TableAllinOnetemp, double& lastEmintemp) {
	using namespace std;
	using namespace arma;
	int ikBTtemp = TableOrderUseOldResultstemp(tableNumnew, 0);
	int id_settemp = TableOrderUseOldResultstemp(tableNumnew, 1);
	//find information for filenames to be loaded
	int tableNumtemp = TableOrdertempInversetemp(ikBTtemp, id_settemp);
	int roundNumtemp = TableAllinOnetemp(tableNumtemp, 14);
	std::string TableRoundNumold = "t" + std::to_string(tableNumtemp) + "r" + std::to_string(roundNumtemp);
	int save_cycle1 = TableAllinOnetemp(tableNumtemp, 15);
	string toloadDensity = prefix + "SeaFillingBG" + TableRoundNumold + to_string((int)save_cycle1) + ".csv";
	mainout << "use density from file: " << toloadDensity << endl;
	if (!preTableSeaDensity.load(toloadDensity, arma::auto_detect)) {
		mainout << "file not found: " << toloadDensity << endl;
	}
	else {
		mainout << "file sucess: " << toloadDensity << endl;
	}
	preTableSeaDensity = fillingToAbs2(preTableSeaDensity);
	string toloadVT = prefix + "VT" + TableRoundNumold + to_string((int)save_cycle1) + ".csv";
	if (!preTableVT.load(toloadVT, arma::auto_detect)) {
		mainout << "file not found: " << toloadVT << endl;
	}
	else {
		mainout << "file sucess: " << toloadVT << endl;
	}
	//information on bounds of eigenenergy
    //(almost when fix kBT, exactly when fix second parameter)the same Hamiltonian will be diagonalized, only the temperature will be lower
	energylowtemp = TableAllinOnetemp(tableNumtemp, 20);
	energyhightemp = TableAllinOnetemp(tableNumtemp, 21);
	lastEmintemp = TableAllinOnetemp(tableNumtemp, 32);//only update, not used here
	double bandwidth = energyhightemp - energylowtemp;
	mainout << "energylowtemp-pre: " << energylowtemp << " energyhightemp-pre: " << energyhightemp << endl;
	double kBTold = MyTemperatureSet(ikBTtemp);
	energylowtemp = energylowtemp - 0.1 * bandwidth - 1.0 * occuMinimumEkBTSecond * kBTold;
	//energyhightemp = energyhightemp + 0.1 * bandwidth + 1.0 * occuMinimumEkBTSecond * kBTold;
	energyhightemp = energyhightemp + 0.02 + 1.0 * occuMinimumEkBTSecond * kBTold;
	/*yayun 20240219, to avoid segmentation fault due to large wavefunction storage; should be fine if table parameters change is small */
	/*
	energylowtemp = energylowtemp;
	energyhightemp = energyhightemp;
	*/
	mainout << "energylowtemp: " << energylowtemp << ", energyhightemp: " << energyhightemp << endl;
	//information on Nsubspace
	/*
	int Nsubspaceold = TableAllinOnetemp(tableNumtemp, 30);
	Nsubspacetemp = Nsubspaceold;
	mainout << "Nsubspace-pre: " << Nsubspaceold << ", Nsubspacetemp: " << Nsubspacetemp << endl;
	*/
	int Neetemp = TableAllinOnetemp(tableNumtemp, 22);//number of eigenstates obtained from diagonalization
	Nsubspacetemp = round((double)((energyhightemp - energylowtemp) / bandwidth) * Neetemp) * 2.0;//estimated number of eigenstates obtained from diagonalization	
	mainout << "Neetemp-pre: " << Neetemp << " Nsubspacetemp: " << Nsubspacetemp << endl;

	return true;
}


//loading old results when starting a new tableNum that has not been started before, only want to make sure I start with previous density and VT, and energy bounds
/*loadPreviousTableResults almost same as loadPreviousRoundResults, apart from the reference parameter name: tableNumnew, and its usage*/
int loadPreviousRoundResults(const int& tableNumtemp, arma::mat& preTableSeaDensity, arma::mat& preTableVT, double& energylowtemp, double& energyhightemp, int& Nsubspacetemp,
	const arma::umat& TableOrdertemp, const arma::mat& TableAllinOnetemp, double& lastEmintemp) {
	using namespace std;
	using namespace arma;
	int ikBTtemp = TableOrdertemp(tableNumtemp, 0);
	int id_settemp = TableOrdertemp(tableNumtemp, 1);
	//find information for filenames to be loaded
	int roundNumtemp = TableAllinOnetemp(tableNumtemp, 14);//check by iterLast in the following if this is really useful
	int iterLast=TableAllinOnetemp(tableNumtemp, 10);

	if (iterLast <= 0) {//this round does not really start
		if (roundNumtemp == 0){
			return 1;//no useful data in this tablenum, 
		}
		else {
			roundNumtemp -= 1;
		}
	}	
	std::string TableRoundNumold = "t" + std::to_string(tableNumtemp) + "r" + std::to_string(roundNumtemp);
	int save_cycle1 = TableAllinOnetemp(tableNumtemp, 15);
	string toloadDensity = prefix + "SeaFillingBG" + TableRoundNumold + to_string((int)save_cycle1) + ".csv";
	mainout << "use density from file: " << toloadDensity << endl;
	if (!preTableSeaDensity.load(toloadDensity, arma::auto_detect)) {
		mainout << "file not found: " << toloadDensity << endl;
	}
	else {
		mainout << "file sucess: " << toloadDensity << endl;
	}
	preTableSeaDensity = fillingToAbs2(preTableSeaDensity);
	string toloadVT = prefix + "VT" + TableRoundNumold + to_string((int)save_cycle1) + ".csv";
	if (!preTableVT.load(toloadVT, arma::auto_detect)) {
		mainout << "file not found: " << toloadVT << endl;
	}
	else {
		mainout << "file sucess: " << toloadVT << endl;
	}
	//information on bounds of eigenenergy
	//(almost when fix kBT, exactly when fix second parameter)the same Hamiltonian will be diagonalized, only the temperature will be lower
	energylowtemp = TableAllinOnetemp(tableNumtemp, 20);
	energyhightemp = TableAllinOnetemp(tableNumtemp, 21);
	lastEmintemp = TableAllinOnetemp(tableNumtemp, 32);//only update, not used here
	double bandwidth = energyhightemp - energylowtemp;
	mainout << "energylowtemp-pre: " << energylowtemp << " energyhightemp-pre: " << energyhightemp << endl;
	double kBTold = MyTemperatureSet(ikBTtemp);
	energylowtemp = energylowtemp - 0.1 * bandwidth - 1.0 * occuMinimumEkBTSecond * kBTold;
	energyhightemp = energyhightemp + 0.1 * bandwidth + 1.0 * occuMinimumEkBTSecond * kBTold;
	/*
	energylowtemp = energylowtemp;
	energyhightemp = energyhightemp;
	*/
	mainout << "energylowtemp: " << energylowtemp << ", energyhightemp: " << energyhightemp << endl;	
	//information on Nsubspace
	/*
	int Nsubspaceold = TableAllinOnetemp(tableNumtemp, 30);
	Nsubspacetemp = Nsubspaceold;
	mainout << "Nsubspace-pre: " << Nsubspaceold << ", Nsubspacetemp: " << Nsubspacetemp << endl;
	*/
	int Neetemp = TableAllinOnetemp(tableNumtemp, 22);//number of eigenstates obtained from diagonalization
	Nsubspacetemp = round((double)((energyhightemp - energylowtemp)/bandwidth)* Neetemp)*2.0;//estimated number of eigenstates obtained from diagonalization	
	mainout << "Neetemp-pre: " << Neetemp << " Nsubspacetemp: " << Nsubspacetemp << endl;
	
	return 0;
}


void assignTableAllinOne(arma::mat& TableAllinOnetemp,const int& tableNumtemp,const double& fermiEnergytemp/*1*/, const double& Etotaltemp/*2*/, const double& KEtemp/*3*/,
	const double& EHartreetemp/*4*/, const double& Eexttemp/*5*/, const double& Excsumtemp/*6*/, const int& itertemp/*7*/, const double& relDensityDifftemp/*8*/, 
	const double& relDensityDiffAimtemp/*9*/, const int& save_cycle1temp/*10*/, const double& Energylowtemp/*11*/, const double& Energyhightemp/*12*/,
	const int& fixNumtemp/*13*/, const int& Nsubspacetemp/*14*/, const double& totalPotentialMinimumtemp/*15*/, const double& lastEmintemp/*16*/,
	const double& nee_wantdecimaltemp/*17*/, const double& Neetemp/*18*/){
		//iterative update tables
		TableAllinOnetemp(tableNumtemp, 0) = fermiEnergytemp; //0FermiEnergy
		//TableAllinOnetemp(tableNumtemp, 1) = 0.0; //1TotalAngualr
		TableAllinOnetemp(tableNumtemp, 2) = Etotaltemp;//2TotalEigenEnergy
		TableAllinOnetemp(tableNumtemp, 3) = KEtemp;//3TotalKinetic
		//TableAllinOnetemp(tableNumtemp, 4) = 0.0;//4TotalKohnShamPotentialEnergy
		//TableAllinOnetemp(tableNumtemp, 5) = 0.0;//5TotalDensityEnergy
		TableAllinOnetemp(tableNumtemp, 6) = EHartreetemp;//6ElectronHartreeEnergy
		TableAllinOnetemp(tableNumtemp, 7) = Eexttemp;//7DiskPotentialEnergy
		//TableAllinOnetemp(tableNumtemp, 8) = 0.0;//8KineticDensityEnergy
		TableAllinOnetemp(tableNumtemp, 9) = Excsumtemp;//9Exchange-correlation-Energy
		TableAllinOnetemp(tableNumtemp, 10) = itertemp;//10iterationNum(this iter also ends here), can be -2;
		TableAllinOnetemp(tableNumtemp, 11) = relDensityDifftemp;//11OverlapRatio
		TableAllinOnetemp(tableNumtemp, 12) = relDensityDiffAimtemp;//12OverlapBetweenAims
		//TableAllinOnetemp(tableNumtemp, 13) = Bigiter;//13RhoBigNum
		//TableAllinOnetemp(tableNumtemp, 14) = roundNum;//14roundNum, only set at the very beginning for-loop of tableNum
		TableAllinOnetemp(tableNumtemp, 15) = save_cycle1temp;//15savecycle1
		//TableAllinOnetemp(tableNumtemp, 16) = save_cycle2;//16savecycle2
		//TableAllinOnetemp(tableNumtemp, 17) = 0.0;//17tableNumStartStage
		//TableAllinOnetemp(tableNumtemp, 18) = iter;//18iterConverge
		//TableAllinOnetemp(tableNumtemp, 19) = -1;//19ConvergeFlag, Converged=1; has been started=0; finished but not converged=-2;default=-1;
		TableAllinOnetemp(tableNumtemp, 20) = Energylowtemp;//20Energylow
		TableAllinOnetemp(tableNumtemp, 21) = Energyhightemp;//21Energyhigh
		//TableAllinOnetemp(tableNumtemp, 22) = Nee;//22number of orbitals  really get
		//TableAllinOnetemp(tableNumtemp, 23) = 0.0;//23temperature
		//TableAllinOnetemp(tableNumtemp, 24) = 0.0;//24magnetic field or MySecondSet
		//TableAllinOnetemp(tableNumtemp, 25) = Lx;//25Lx
		//TableAllinOnetemp(tableNumtemp, 26) = Ly;//26Ly
		//TableAllinOnetemp(tableNumtemp, 27) = Nrx;//27Nrx
		//TableAllinOnetemp(tableNumtemp, 28) = Nry;//28Nry
		TableAllinOnetemp(tableNumtemp, 29) = fixNumtemp;//29Nry
		TableAllinOnetemp(tableNumtemp, 30) = Nsubspacetemp;//30Nsubspace
		TableAllinOnetemp(tableNumtemp, 31) = totalPotentialMinimumtemp;//31totalPotentialMinimum,averaged from the lowest 20?
		TableAllinOnetemp(tableNumtemp, 32) = lastEmintemp;//32 lowest eigen energy of the obtanied orbitals
		TableAllinOnetemp(tableNumtemp, 33) = nee_wantdecimaltemp;//33 total particle number
		TableAllinOnetemp(tableNumtemp, 34) = Neetemp;//22number of orbitals  really used in occu etc
		/*
		//iterative update tables, put it here in case useful for the future
		TableAllinOne(tableNum, 0) = fermiEnergy; //0FermiEnergy
		//TableAllinOne(tableNum, 1) = 0.0; //1TotalAngualr
		TableAllinOne(tableNum, 2) = Etotal;//2TotalEigenEnergy
		TableAllinOne(tableNum, 3) = KE;//3TotalKinetic
		//TableAllinOne(tableNum, 4) = 0.0;//4TotalKohnShamPotentialEnergy
		//TableAllinOne(tableNum, 5) = 0.0;//5TotalDensityEnergy
		TableAllinOne(tableNum, 6) = EHartree;//6ElectronHartreeEnergy
		TableAllinOne(tableNum, 7) = Eext;//7DiskPotentialEnergy
		//TableAllinOne(tableNum, 8) = 0.0;//8KineticDensityEnergy
		TableAllinOne(tableNum, 9) = Excsum;//9Exchange-correlation-Energy
		TableAllinOne(tableNum, 10) = iter;//10iterationNum(this iter also ends here), can be -2;
		TableAllinOne(tableNum, 11) = relDensityDiff;//11OverlapRatio
		TableAllinOne(tableNum, 12) = relDensityDiffAim;//12OverlapBetweenAims
		//TableAllinOne(tableNum, 13) = Bigiter;//13RhoBigNum
		//TableAllinOne(tableNum, 14) = roundNum;//14roundNum, only set at the very beginning for-loop of tablenum
		TableAllinOne(tableNum, 15) = save_cycle1;//15savecycle1
		//TableAllinOne(tableNum, 16) = save_cycle2;//16savecycle2
		//TableAllinOne(tableNum, 17) = 0.0;//17tableNumStartStage
		//TableAllinOne(tableNum, 18) = iter;//18iterConverge
		//TableAllinOne(tableNum, 19) = -1;//19ConvergeFlag, Converged=1; has been started=0; finished but not converged=-2;default=-1;
		TableAllinOne(tableNum, 20) = Energylow;//20Energylow
		TableAllinOne(tableNum, 21) = Energyhigh;//21Energyhigh
		//TableAllinOne(tableNum, 22) = Nee;//22number of orbitals
		//TableAllinOne(tableNum, 23) = 0.0;//23temperature
		//TableAllinOne(tableNum, 24) = 0.0;//24magnetic field or MySecondSet
		//TableAllinOne(tableNum, 25) = Lx;//25Lx
		//TableAllinOne(tableNum, 26) = Ly;//26Ly
		//TableAllinOne(tableNum, 27) = Nrx;//27Nrx
		//TableAllinOne(tableNum, 28) = Nry;//28Nry
		TableAllinOne(tableNum, 29) = fixNum;//29Nry
		TableAllinOne(tableNum, 30) = Nsubspace;//30Nsubspace
		*/
}


void saveHamToFile(const string& prefixtemp, const string& TableRoundNumtemp, const int& save_cycle1temp) {
	using namespace arma;
	//INFO(Nx);
	//INFO(Ny);
	//INFO(Nxy);
	//INFO(ham.nnz);
	{
		ofstream f(prefix + "matHam" + TableRoundNumtemp + to_string((int)save_cycle1temp), std::ofstream::trunc);
		f << std::setprecision(18);
		for (auto i = 0; i < ham.nnz; i++) {
			f << ham.rows1[i] << '\t' << ham.cols1[i] << '\t' << real(ham.vals[i]) << '\t'
			  << imag(ham.vals[i]) << endl;
		}
		f.close();
	}
	{
		ofstream f(prefix + "matHam" + TableRoundNumtemp + to_string((int)save_cycle1temp) + ".ana", std::ofstream::trunc);
		f << std::setprecision(18);
		f << "#row\tcol\tv.real\tv.imag\tix\tiy\tjx\tjy" << endl;// This is the header line or format of writtings
		for (auto i = 0; i < ham.nnz; i++) {
			int ix, iy, jx, jy;
			gridFromIndex(ham.rows[i], ix, iy);
			gridFromIndex(ham.cols[i], jx, jy);
			f << ham.rows1[i] << '\t' << ham.cols1[i] << '\t' << real(ham.vals[i]) << '\t'
			  << imag(ham.vals[i]) << "\t" << ix << '\t' << iy << '\t' << jx << '\t'<<jy<<endl;
		}
		f.close();
	}
}



void saveDataToFile(const string& prefixtemp, const string& TableRoundNumtemp, const int& save_cycle1temp,const arma::mat & SeaDensitytemp, 
	const arma::mat& newSeaDensitytemp,
	const arma::mat& oldVTtemp, const arma::mat& A_x_mattemp, const arma::mat& A_y_mattemp, const arma::mat& VCoulombtemp, const arma::mat totalPotentialVxctemp
	) {
	abs2ToFilling(SeaDensitytemp).save(prefixtemp + "SeaFillingBG" + TableRoundNumtemp + to_string((int)save_cycle1temp) + ".csv", arma::csv_ascii);
	abs2ToFilling(newSeaDensitytemp).save(prefixtemp + "SeaFillingBGnew" + TableRoundNumtemp + to_string((int)save_cycle1temp) + ".csv", arma::csv_ascii);
	oldVTtemp.save(prefixtemp + "VT" + TableRoundNumtemp + to_string((int)save_cycle1temp) + ".csv", arma::csv_ascii);
	A_x_mattemp.save(prefixtemp + "Ax" + TableRoundNumtemp + to_string((int)save_cycle1temp) + ".csv", arma::csv_ascii);
	A_y_mattemp.save(prefixtemp + "Ay" + TableRoundNumtemp + to_string((int)save_cycle1temp) + ".csv", arma::csv_ascii);
	VCoulombtemp.save(prefixtemp + "VCoulomb" + TableRoundNumtemp + to_string((int)save_cycle1temp) + ".csv", arma::csv_ascii);
	totalPotentialVxctemp.save(prefixtemp + "gradVxc" + TableRoundNumtemp + to_string((int)save_cycle1temp) + ".csv", arma::csv_ascii);
	/*
			abs2ToFilling(SeaDensity).save(prefix + "SeaFillingBG" + TableRoundNum + to_string((int)save_cycle1) + ".csv", arma::csv_ascii);
			abs2ToFilling(newSeaDensity).save(prefix + "SeaFillingBGnew" + TableRoundNum + to_string((int)save_cycle1) + ".csv", arma::csv_ascii);
			oldVT.save(prefix + "VT" + TableRoundNum + to_string((int)save_cycle1) + ".csv", arma::csv_ascii);
			A_x_mat.save(prefix + "Ax" + TableRoundNum + to_string((int)save_cycle1) + ".csv", arma::csv_ascii);
			A_y_mat.save(prefix + "Ay" + TableRoundNum + to_string((int)save_cycle1) + ".csv", arma::csv_ascii);
			VCoulomb.save(prefix + "VCoulomb" + TableRoundNum + to_string((int)save_cycle1) + ".csv", arma::csv_ascii);
			totalPotentialVxc.save(prefix + "gradVxc" + TableRoundNum + to_string((int)save_cycle1) + ".csv", arma::csv_ascii);
	*/
}

/* to make a table for reference of annealling, TableOrder controls the order of parameter usage, 
TableOrderUseOldResults controls what old results to be used as initial input; 0-based*/
void MakeTableOrderDimension(int const TemperatureNum, int const MySecondSetNum, arma::umat& TableOrdertemp, arma::umat& TableOrdertempInversetemp, 
	arma::umat& TableOrderUseOldResultstemp){
	using namespace std;
	using namespace arma;
	int OrderNum;
	OrderNum = TemperatureNum * MySecondSetNum;
	TableOrdertemp.zeros(OrderNum, 2);
	TableOrderUseOldResultstemp.zeros(OrderNum, 2);
	TableOrdertempInversetemp.zeros(TemperatureNum, MySecondSetNum);
	OrderNum = 0;
	for (int i = 0; i <= TemperatureNum - 1; i++) {
		for (int j = 0; j <= MySecondSetNum - 1; j++) {
			TableOrdertemp(OrderNum, 0) = i;
			TableOrdertemp(OrderNum, 1) = j;//always go along MySecondSetNum for a fixed temperature, TableOrder records the i and j indices of temperature and secondset used in OrderNum^th step; 
			if (i == 0) {
				TableOrderUseOldResultstemp(OrderNum, 0) = i;
				if(j == 0){
					TableOrderUseOldResultstemp(OrderNum, 1) = 0;// very first step is special
				}
				else
				{
					TableOrderUseOldResultstemp(OrderNum, 1) = j - 1;// at the first temperature, use previous MySecondSetNum as initial input
				}
			}
			else if(i > 0){
				TableOrderUseOldResultstemp(OrderNum, 0) = i - 1;// at later temperatures, use the previous temperature but the same MySecondSetNum as initial input
				TableOrderUseOldResultstemp(OrderNum, 1) = j;
			}
			TableOrdertempInversetemp(i, j) = OrderNum;
			OrderNum++;
		}
	}
}


/*ConvergenceTest, the input arrays effectively start from 1*/ //NOT tested
//ConvergenceTest(iter, Bigiter, AllbigDiff, AlldensityDiff, Allnewparticle)
bool ConvergenceTest(int& nrun, int& Bigiter, arma::vec& AllbigDiff, arma::vec& AlldensityDiff, arma::vec& Allnewparticle) {/*arma::vec AlltotalEnergy*/
	//lookBackStepNum,Bigiter is the index number of the first effective step and Bigiter+1 is the corresponding number of iter
	using namespace std;
	using namespace arma;
	if (nrun > lookBackStepNum && nrun-lookBackStepNum > Bigiter && Bigiter!= 0) {
		arma::vec cutAllbigDiff=AllbigDiff.cols(nrun- lookBackStepNum, nrun - 1);
		arma::vec cutAlldensityDiff= AlldensityDiff.cols(nrun- lookBackStepNum, nrun - 1);
		arma::vec cutAllnewparticle = Allnewparticle.cols(nrun- lookBackStepNum, nrun - 1);
		if (cutAllbigDiff.max() - cutAllbigDiff.min() < BigNumThreshold && cutAlldensityDiff.max()< relDensityThreshold && 
			cutAllnewparticle.max()- cutAllnewparticle.min()< chemicalNumThreshold){
			//to add monotonicity comparison in the future
			return true;
		}
		else {
			return false;
		}
	}else{//not enough steps
		return false;
	}
}


/*ConvergenceRecent, the input arrays effectively start from 1*/
//
bool ConvergenceRecent(int& nrun, int& Bigiter, double& KE, arma::vec& cutAllbigDiff, arma::vec& cutAlldensityDiff, arma::vec& cutAllVTDiff, arma::vec& cutAllkineticEnergy,
	arma::vec& cutAllnewparticle) {
	//lookBackStepNum,Bigiter is the index number of the first effective step and Bigiter+1 is the corresponding number of iter
	using namespace std;
	using namespace arma;
	if (nrun > lookBackStepNum && nrun - lookBackStepNum > Bigiter && Bigiter != 0) {
		if (cutAllbigDiff.max() - cutAllbigDiff.min() < BigNumThreshold /*this also includes the particle change*/ 
			&& cutAlldensityDiff.max() < relDensityThreshold && cutAllVTDiff.max() < relDensityThreshold
			&& abs(cutAllnewparticle.max() - cutAllnewparticle.min()) < chemicalNumThreshold/*particle number fluctuation magnitude*/
			&& abs(cutAllkineticEnergy.max() - cutAllkineticEnergy.min())*0.0 < abs(KineticThreshold*KE)/*particle number fluctuation magnitude*/) {
			//to add monotonicity comparison in the future
			return true;
		}
		else {
			return false;
		}
	}
	else {//not enough steps
		return false;
	}
}


/*ConvergenceRecent, the input arrays effectively start from 1*/
//
bool ConvergenceRecentNoBig(int& nrun, int& Bigiter, double& KE, arma::vec& cutAllbigDiff, arma::vec& cutAlldensityDiff, arma::vec& cutAllVTDiff, arma::vec& cutAllkineticEnergy,
	arma::vec& cutAllnewparticle) {
	//lookBackStepNum,Bigiter is the index number of the first effective step and Bigiter+1 is the corresponding number of iter
	using namespace std;
	using namespace arma;
	if (nrun > lookBackStepNum) {
		if (cutAllbigDiff.max() - cutAllbigDiff.min() < BigNumThreshold /*this also includes the particle change*/
			&& cutAlldensityDiff.max() < relDensityThreshold && cutAllVTDiff.max() < relDensityThreshold
			&& abs(cutAllnewparticle.max() - cutAllnewparticle.min()) < chemicalNumThreshold/*particle number fluctuation magnitude*/
			//&& abs(cutAllkineticEnergy.max() - cutAllkineticEnergy.min()) * 0.0 < abs(KineticThreshold * KE)/*particle number fluctuation magnitude*/
			) {
			//to add monotonicity comparison in the future
			return true;
		}
		else {
			return false;
		}
	}
	else {//not enough steps
		return false;
	}
}



/*find potential minimum*/
double FindPotentialMinimum(const arma::mat& totalPotentialtemp, const int& averNum) {
	arma::vec totalPotentialcolumn = sort(totalPotentialtemp.as_col());
	return mean(totalPotentialcolumn.rows(1, averNum));
}
/*check file accessibility*/
inline bool exists_test(const std::string& name) {
	if (FILE* file = fopen(name.c_str(), "r")) {
		fclose(file);
		return true;
	}
	else {
		return false;
	}
}


/*TODO: yayun 231231. find eigenstate of kinetic plus barrier potential within the [0, 0.11*3] energy window, 
then find thier expection in others potential terms*/
double FindKinetic(const arma::mat& totalPotentialtemp, const int& averNum) {
	arma::vec totalPotentialcolumn = sort(totalPotentialtemp.as_col());
	return mean(totalPotentialcolumn.rows(1, averNum));
}



/* Record occupations of Table parameters in binary. */
void recordoccuTablebinary(const string& filename, const int& nrun, const int& nee) {
	using namespace std;
	ofstream file(filename, nrun == 0 ? (ios_base::out | ios::binary) : (ios_base::out | ios_base::app | ios::binary));
	auto ptr = occupationNumber;
	size_t len = sizeof(ptr[0]) * nee;
	file.write((char*)ptr, len);
	file.close();
}

/* Record ees of Table parameters in binary. */
void recordEigenTablebinary(const string& filename, const int& nrun, const int& nee) {
	using namespace std;
	ofstream file(filename, nrun == 0 ? (ios_base::out | ios::binary) : (ios_base::out | ios_base::app | ios::binary));
	auto ptr = Eigenvals;
	size_t len = sizeof(ptr[0]) * nee;
	file.write((char*)ptr, len);
	file.close();
}


/*yayun 240130. sort an array*/
void SortPointerArrayIntoArmaArray(arma::rowvec& sorted_Eigen, int nee_get, const double* ees) {
	/*
	for (int it = nee_get-1; it > nee_get-20; --it) { 
		cout << ees[it] << "; "; 
	}	
	cout << endl;
	*/
	for (int it = 0; it < nee_get; ++it) { 
		sorted_Eigen[it] = ees[it]; 
	}
	sorted_Eigen = sort(sorted_Eigen, "ascend", 1);
	/*
	for (int it = nee_get-1; it > nee_get-20; --it) { 
		cout << sorted_Eigen[it] << "; "; 
	}
	cout << endl;
	cout << "sorted_Eigen[0]: " << sorted_Eigen[0] << "sorted_Eigen[nee_get-1]: " << sorted_Eigen[nee_get-1] << endl;
	*/
	return;
}












#endif
