#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

extern int build_bath(int n, gsl_complex ***Green_bath, PNSpara egv);
extern int myrank ;
extern gsl_complex ***Tbath ;
extern gsl_complex ****Greennew ;
extern PNSpara fegv;
using namespace std;

int tNB=2*NB;

class infoHybblock { 
	public : 
		int nblock , nbasis , ispin ;
		vector< vector<int> > blockIndArr ; 
		vector< vector<int> > bathBlockIndArr ; 
		vector<int> nbathBlockArr ; 
		vector<int> blockSymArr ; 
		vector<int> blockSymIndArr ; 

		infoHybblock(void) {}
		infoHybblock( char *fname ) {
			ifstream fin ;
			fin.open( fname );
			while( !fin.eof() ) {
				char line[1024] ;
				fin.getline( line, 1024 );
				string str ;
				stringstream ss ;
				ss << line ;
				//cout << line <<endl;

				str = ss.str() ;
				if( str.size()>0 ) {
					string numberstr(str) ; numberstr.erase( numberstr.find(":") , numberstr.size() );
					int dum=-1;
					//cout << "Found at  "<< str.find(":")<<" / "<<str.size()<<"\n";
					dum = atoi( numberstr.c_str() );
					//cout << "\tAFTER ERASE  : "<< numberstr <<endl;
					//cout << "\tDUM  : "<< dum <<endl;
					//cout << "\tTEST : "
					//	<<";"<< str.find("total nbasis in a impurity")
					//	<<";"<< str.find("number of blocks within block-diagonal")
					//	<<";"<< str.find("turn-on (1) / turn-off (0) for spin-block")
					//	<<";"<< str.find("n-th block, its indices") << endl;

					if( str.find("total nbasis in a impurity")<string::npos )				nbasis = dum ;
					else if( str.find("number of blocks within block-diagonal")<string::npos ) {		nblock = dum ; blockIndArr.resize( nblock ); } //cout <<"\tNBLOCK : "<<nblock<<endl;}
					else if( str.find("turn-on (1) / turn-off (0) for spin-block")<string::npos )		ispin  = dum ;
					else if( str.find("n-th block, its indices")<string::npos ) {
						string w ;
						vector<int> word ;
						for(stringstream sts(numberstr);(sts>>w);) {
							//cout << dum << " W : " << w << endl ;
							word.push_back( atoi(w.c_str()) );
						}
						blockIndArr[dum].resize( word.size()-1 );
						word.erase( word.begin() );
						blockIndArr[dum].swap( word ) ;
					}
					else if( str.find("number of bath-orbital in each block")<string::npos ) {
						string w ;
						vector<int> word ;
						for(stringstream sts(numberstr);(sts>>w);) {
							//cout << dum << " W : " << w << endl ;
							word.push_back( atoi(w.c_str()) );
						}
						nbathBlockArr.resize( word.size() );
						nbathBlockArr.swap( word );
						if( abs( nbathBlockArr.size()- nblock ) )	cout<<"Warnning :: size of nbathBlockArr does not match with nblock.\n";

						bathBlockIndArr.resize( nblock );
						int itot = 0 ;
						for( int i=0 ; i<nblock ; i++ ){
							for( int j=0 ; j<nbathBlockArr[i] ; j++ ){
								bathBlockIndArr[i].push_back( j + itot );
							}
							itot += nbathBlockArr[i] ;
						}
					}
					else if( str.find("block indices to be symmetrized")<string::npos ) {
						string w ;
						vector<int> word ;
						for(stringstream sts(numberstr);(sts>>w);) {
							//cout << dum << " W : " << w << endl ;
							word.push_back( atoi(w.c_str()) );
						}
						blockSymArr.resize( word.size() );
						blockSymArr.swap( word );
					}
					else 											cout << "Unknown variables. ; "<< str << endl;
				}
			}
			blockSymIndArr.resize(nblock);
			for( unsigned int i=0 ; i<blockSymArr.size() ; i++ ){
				blockSymIndArr[ blockSymArr[i] ] = 1 ;
			}
		}
		void show() {
			cout<< "-------------------------\n";
			cout<<"\tnblock , nbasis , ispin ; "<< nblock <<" "<<  nbasis <<" "<<  ispin << endl ;
			cout<< "blockIndArr.size : " << blockIndArr.size() << endl;
			for( unsigned int i=0 ; i<blockIndArr.size() ; i++ ) {
				cout<<"\t"<<i<<"-th block orb's : ";
				for( unsigned int j=0 ; j<blockIndArr[i].size() ; j++ ) 
					cout << blockIndArr[i][j] <<" ";
				cout << endl ;
			}
			cout<< "nbathBlockArr : " ;
			for( unsigned int i=0 ; i<nbathBlockArr.size() ; i++ ){	cout << nbathBlockArr[i] <<" "; } cout << endl ;
			cout<< "bathBlockIndArr.size : " << bathBlockIndArr.size() << endl;
			for( unsigned int i=0 ; i<bathBlockIndArr.size() ; i++ ) {
				cout<<"\t"<<i<<"-th block bath's : ";
				for( unsigned int j=0 ; j<bathBlockIndArr[i].size() ; j++ ) 
					cout << bathBlockIndArr[i][j] <<" ";
				cout << endl ;
			}

			cout<< "blockSymArr.size : " << blockSymArr.size() << endl;
			for( unsigned int i=0 ; i<blockSymArr.size() ; i++ ) {
				cout << "\t"<< blockSymArr[i] << "-th block will be symmetrized." <<endl ;
			}
			cout<< "blockSymIndArr : " ;
			for( unsigned int i=0 ; i<blockSymIndArr.size() ; i++ ) {
				cout << blockSymIndArr[i] << " ";
			}
			cout <<endl;

			cout<< "-------------------------\n";
			//vector<int> vint( blockIndArr , blockIndArr+nblock ); 
			//int iv=0;
			//for(vector<int>::iterator it = vint.begin() ; it!=vint.end() ; ++it) { cout<<iv<<"-th block : "<< *it << endl ; iv++; }
		}
		void read( char *fname ) {
			ifstream fin ;
			fin.open( fname );
			while( !fin.eof() ) {
				string s, w;
				getline( fin, s );
				vector<string> word;
				cout << "W : \n";
				for(stringstream sts(s);(sts>>w);){
					word.push_back(w);
					cout << w << " " ;
				}
				cout << "W.end \n";
				for(vector<string>::iterator it = word.begin();it!=word.end();++it) { 
					if( strcmp((*it).c_str(),":")==0 ) {
						break ;
					}
					cout<<*it<< " " ;
				}
				cout << endl;
			}
		}
		int getNp(int block){
			if( blockSymIndArr[block]<1 )
				return 2*nbathBlockArr[block]*blockIndArr[block].size() + nbathBlockArr[block];
			else
                                return 2*nbathBlockArr[block]/blockIndArr[block].size() + nbathBlockArr[block];
		}
};

int convert_egv_to_p_block(PNSpara egv, double p[], int sp, int block, int nInd, int *indArr, int nIndBath, int *indBathArr ){
	for(int i=0; i<nIndBath; i++)		p[i] = egv.egbath[indBathArr[i]];
	for(int i=0; i<nInd; i++) {
		int offset = nIndBath + nIndBath*i;
		int mu = indArr[i];
		for(int j=0; j<nIndBath; j++){
			int l  = indBathArr[j];
			p[offset + j ]			= GSL_REAL( egv.hybrid[mu][l] );
			p[offset + j + nIndBath*nInd]	= GSL_IMAG( egv.hybrid[mu][l] );
		}
	}
	return 0;
}
int convert_p_to_egv_block(PNSpara egv, double p[], int sp, int block, int nInd, int *indArr, int nIndBath, int *indBathArr ){
	for(int i=0; i<nIndBath; i++)		egv.egbath[indBathArr[i]] = p[i];
	for(int i=0; i<nInd; i++) { 
		int offset = nIndBath + nIndBath*i;
		int mu = indArr[i];
		for(int j=0; j<nIndBath; j++){
			int l  = indBathArr[j];
			egv.hybrid[mu][l] =
				gsl_complex_rect(
						p[offset + j],
						p[offset + j + nIndBath*nInd] 
						);
		}
	}
	return 0;
}
int convert_egv_to_p_blockSym(PNSpara egv, double p[], int sp, int block, int nInd, int *indArr, int nIndBath, int *indBathArr){
	int nIndBathEnergy = nIndBath/nInd;
	for(int i=0; i<nIndBathEnergy; i++)		p[i] = egv.egbath[ indBathArr[i*nInd] ];
	//for(int i=0; i<nIndBathEnergy; i++)		cout << "p_el "<< indBathArr[i*nInd] << " : "<< p[i]<<endl;
	for(int i=0; i<nIndBathEnergy; i++) {
		int mu = indArr[0];
		int l  = indBathArr[ i*nInd ];
		p[ nIndBathEnergy + i ]			= gsl_complex_abs( egv.hybrid[mu][l] );
		//cout << "p_v "<< mu<<" "<<l<<" : "<< p[ nIndBathEnergy + i ] <<endl ;
	}

	int nIndBathEnergy2 = 2*nIndBathEnergy;
	for(int i=0; i<nIndBathEnergy; i++){
		int offset2	= nIndBathEnergy2 + i*nInd;
		for(int imu=0; imu<nInd; imu++) {
			int mu	= indArr[ imu ];
			int l	= indBathArr[ i*nInd + imu ];
			p[ offset2 + imu ]		= gsl_complex_arg( egv.hybrid[mu][l] );
		}
	}
	//cout<< "PtoEGV-------------------------------\n";
	//for(int l=0; l<tNB; l++){ 
	//	printf("%12.9f\t", egv.egbath[ l ]);
	//	for(int mu=0; mu<tNC; mu++){
	//		printf("%12.9f\t", GSL_REAL(egv.hybrid[mu][l]) );
	//		printf("%12.9f\t", GSL_IMAG(egv.hybrid[mu][l]) );
	//	}
	//	cout << endl;
	//}
	//cout<< "-------------------------------------\n";
	return 0;
}
int convert_p_to_egv_blockSym(PNSpara egv, double p[], int sp, int block, int nInd, int *indArr, int nIndBath, int *indBathArr){
	int nIndBathEnergy = nIndBath/nInd;

	for(int i=0; i<nIndBathEnergy; i++)
		for(int imu=0; imu<nInd; imu++) {
			egv.egbath[ indBathArr[i*nInd+imu] ]	= p[i] ;
			//cout << "p_el "<< indBathArr[i*nInd+imu] << " : "<< p[i]<<endl;
		}

	//cout << "nIndBath : "<< nIndBath<<endl;
	int nIndBathEnergy2 = 2*nIndBathEnergy;
	for(int i=0; i<nIndBathEnergy; i++) {
		int offset	= nIndBathEnergy;
		int offset2	= nIndBathEnergy2 + i*nInd;
		double Vabs	= p[offset +i];
		for(int imu=0; imu<nInd; imu++) {
			//cout << "i, imu, i*nInd+imu, indBath  : "<< i<<", "<<imu<<", "<< i*nInd + imu <<", "<< indBathArr[i*nInd+imu]<<endl;
			int mu		= indArr[ imu ],
			    l		= indBathArr[ i*nInd + imu ];
			double Varg	= p[ offset2 + imu ];
			egv.hybrid[mu][l]		= gsl_complex_polar( Vabs, Varg );
		}
	}
	//cout<< "EGVtoP-------------------------------\n";
	//for(int l=0; l<tNB; l++){ 
	//	printf("%12.9f\t", egv.egbath[ l ]);
	//	for(int mu=0; mu<tNC; mu++){
	//		printf("%12.9f\t", GSL_REAL(egv.hybrid[mu][l]) );
	//		printf("%12.9f\t", GSL_IMAG(egv.hybrid[mu][l]) );
	//	}
	//	cout << endl;
	//}
	//cout<< "-------------------------------------\n";
	return 0;
}

double my_f_block(double p[], int n, int sp, int block, int nInd, int *indArr, int nIndBath, int *indBathArr ){
	double chisquare = 0, ctemp=0;

	convert_p_to_egv_block(fegv, p, sp, block, nInd, indArr, nIndBath, indBathArr);

	build_bath(n, Tbath, fegv);

#pragma omp parallel for ordered private(ctemp) shared(Greennew, Tbath, Nmax, n, nctot, block, Wwn) default(shared) reduction(+: chisquare)
	for(int i=0; i<Nmax; i++){
		ctemp = 0;
		for(int imu=0; imu<nInd; imu++) {
			int mu = indArr[imu];
			for(int inu=0; inu<nInd; inu++){
				int nu = indArr[inu];
				ctemp += (GSL_REAL(Greennew[n][i][mu][nu])-GSL_REAL(Tbath[i][mu][nu]))*(GSL_REAL(Greennew[n][i][mu][nu])-GSL_REAL(Tbath[i][mu][nu])) +
					(GSL_IMAG(Greennew[n][i][mu][nu])-GSL_IMAG(Tbath[i][mu][nu]))*(GSL_IMAG(Greennew[n][i][mu][nu])-GSL_IMAG(Tbath[i][mu][nu])) ;
			}
		}
#pragma omp ordered
		chisquare += ctemp * WWN;
	}
	return chisquare / Nmax / NC/NC;
}
inline double my_f_block_doublepar(double p[], void *params){
	return my_f_block( p, ((int **)params)[3][0], ((int **)params)[4][0], ((int **)params)[5][0], ((int **)params)[6][0], ((int **)params)[7], ((int **)params)[8][0], ((int **)params)[9] );
}
double my_f_block_gslpar(const gsl_vector *v, void *params){
	int NpBlock = ((int **)params)[10][0];
	double p[NpBlockCap];
	for(int i=0; i<NpBlock; i++)	p[i] = gsl_vector_get(v, i);
	return my_f_block_doublepar( p, params );
}
void my_df_block(const gsl_vector *v, void *params, gsl_vector *df){
	int i,j;

	int **par = (int **)params;
	int *offset	= par[0];
	int *number	= par[1];
	int mystart	= par[2][0], myend = par[2][1];
	int n		= par[3][0];
	int sp		= par[4][0];
	int block	= par[5][0];
	int nInd	= par[6][0];
	int *indArr	= par[7];
	int nIndBath	= par[8][0];
	int *indBathArr	= par[9];
	int NpBlock	= par[10][0];

	double p2plus[NpBlockCap], pplus[NpBlockCap], pminus[NpBlockCap], p2minus[NpBlockCap], h=diff/2.0L;
	double p[NpBlockCap], dp[NpBlockCap];
	for(i=0; i<NpBlock; i++)	p[i] = gsl_vector_get(v, i);

	double h2 = h*2.;
	for(i=mystart; i<myend; i++){
		for(j=0; j<NpBlock; j++){
			p2plus[j]	= p[j];
			pplus[j]	= p[j];
			pminus[j]	= p[j];
			p2minus[j]	= p[j];
		}
		p2plus[i]	+= h2;
		pplus[i]	+= h;
		pminus[i]	-= h;
		p2minus[i]	-= h2;
		dp[i] = (
				-   my_f_block(p2plus,	n, sp, block, nInd, indArr, nIndBath, indBathArr)
				+ 8*my_f_block(pplus,	n, sp, block, nInd, indArr, nIndBath, indBathArr)
				- 8*my_f_block(pminus,	n, sp, block, nInd, indArr, nIndBath, indBathArr) 
				+   my_f_block(p2minus,	n, sp, block, nInd, indArr, nIndBath, indBathArr)
			) / (12.*h);
	}
	MPI_Allgatherv(MPI_IN_PLACE, myend-mystart, MPI_DOUBLE, dp, number, offset, MPI_DOUBLE, MPI_COMM_WORLD);

	for(i=0; i<NpBlock; i++){
		gsl_vector_set( df,  i, dp[i] );
	}
	//printf("myrank = %d, [%d:%d] of %d\n", myrank, mystart, myend, NpBlock);
}
void my_fdf_block(const gsl_vector *x, void *params, double *f, gsl_vector *df){
	*f = my_f_block_gslpar(x, params);
	my_df_block(x, params, df);
}

double my_f_blockSym(double p[], int n, int sp, int block, int nInd, int *indArr, int nIndBath, int *indBathArr ){
	double chisquare = 0, ctemp=0;

	convert_p_to_egv_blockSym(fegv, p, sp, block, nInd, indArr, nIndBath, indBathArr);

	build_bath(n, Tbath, fegv);

#pragma omp parallel for ordered private(ctemp) shared(Greennew, Tbath, Nmax, n, nctot, block, Wwn) default(shared) reduction(+: chisquare)
	for(int i=0; i<Nmax; i++){
		ctemp = 0;
		for(int imu=0; imu<nInd; imu++) {
			int mu = indArr[imu];
			for(int inu=0; inu<nInd; inu++){
				int nu = indArr[inu];
				ctemp += (GSL_REAL(Greennew[n][i][mu][nu])-GSL_REAL(Tbath[i][mu][nu]))*(GSL_REAL(Greennew[n][i][mu][nu])-GSL_REAL(Tbath[i][mu][nu])) +
					(GSL_IMAG(Greennew[n][i][mu][nu])-GSL_IMAG(Tbath[i][mu][nu]))*(GSL_IMAG(Greennew[n][i][mu][nu])-GSL_IMAG(Tbath[i][mu][nu])) ;
			}
		}
#pragma omp ordered
		chisquare += ctemp * WWN;
	}
	return chisquare / Nmax / NC/NC;
}
inline double my_f_block_doubleparSym(double p[], void *params){
	return my_f_blockSym( p, ((int **)params)[3][0], ((int **)params)[4][0], ((int **)params)[5][0], ((int **)params)[6][0], ((int **)params)[7], ((int **)params)[8][0], ((int **)params)[9] );
}
double my_f_block_gslparSym(const gsl_vector *v, void *params){
	int NpBlock = ((int **)params)[10][0];
	double p[NpBlockCap];
	for(int i=0; i<NpBlock; i++)	p[i] = gsl_vector_get(v, i);
	return my_f_block_doubleparSym( p, params );
}
void my_df_blockSym(const gsl_vector *v, void *params, gsl_vector *df){
	int i,j;

	int **par = (int **)params;
	int *offset	= par[0];
	int *number	= par[1];
	int mystart	= par[2][0], myend = par[2][1];
	int n		= par[3][0];
	int sp		= par[4][0];
	int block	= par[5][0];
	int nInd	= par[6][0];
	int *indArr	= par[7];
	int nIndBath	= par[8][0];
	int *indBathArr	= par[9];
	int NpBlock	= par[10][0];

	double p2plus[NpBlockCap], pplus[NpBlockCap], pminus[NpBlockCap], p2minus[NpBlockCap], h=diff/2.0L;
	double p[NpBlockCap], dp[NpBlockCap];
	for(i=0; i<NpBlock; i++)	p[i] = gsl_vector_get(v, i);

	double h2 = h*2.;
	for(i=mystart; i<myend; i++){
		for(j=0; j<NpBlock; j++){
			p2plus[j]	= p[j];
			pplus[j]	= p[j];
			pminus[j]	= p[j];
			p2minus[j]	= p[j];
		}
		p2plus[i]	+= h2;
		pplus[i]	+= h;
		pminus[i]	-= h;
		p2minus[i]	-= h2;
		dp[i] = (
				-   my_f_blockSym(p2plus,	n, sp, block, nInd, indArr, nIndBath, indBathArr)
				+ 8*my_f_blockSym(pplus,	n, sp, block, nInd, indArr, nIndBath, indBathArr)
				- 8*my_f_blockSym(pminus,	n, sp, block, nInd, indArr, nIndBath, indBathArr) 
				+   my_f_blockSym(p2minus,	n, sp, block, nInd, indArr, nIndBath, indBathArr)
			) / (12.*h);
	}
	MPI_Allgatherv(MPI_IN_PLACE, myend-mystart, MPI_DOUBLE, dp, number, offset, MPI_DOUBLE, MPI_COMM_WORLD);

	for(i=0; i<NpBlock; i++){
		gsl_vector_set( df,  i, dp[i] );
	}
	//printf("myrank = %d, [%d:%d] of %d\n", myrank, mystart, myend, NpBlock);
}
void my_fdf_blockSym(const gsl_vector *x, void *params, double *f, gsl_vector *df){
	*f = my_f_block_gslparSym(x, params);
	my_df_blockSym(x, params, df);
}

void minimize_blockAll(int count, PNSpara egv, double *minimal, int *iter, double Converge, int n, int sp, int block){
	infoHybblock iHyb( "inputs/HYBBLOCK" );
	int NpBlock	= iHyb.getNp( block ) ;			if(myrank<1) cout<< "NpBlock "<<block<<" : "<< NpBlock << endl ;
	int nInd	= iHyb.blockIndArr[block].size();
	int indArr[tNC];
	for(int i=0; i<nInd; i++ ) indArr[i] = iHyb.blockIndArr[block][i] + sp ;

	int nIndBath	= iHyb.nbathBlockArr[block];
	int indBathArr[Np];
	for(int i=0; i<nIndBath; i++ ) indBathArr[i] = iHyb.bathBlockIndArr[block][i]*(iHyb.ispin+1) + sp ;
	if(myrank<1){
		printf("\tindices of bath-orbital : ");
	       	for(int i=0; i<nIndBath; i++) printf("%d ",indBathArr[i]);
		printf("\n") ;
	}

	if( NpBlockCap < NpBlock ){ if(myrank<1) printf("Error :: NpBlockCap<NpBlock. %d<%d\n",NpBlockCap,NpBlock); exit(1); }
	double p[NpBlockCap];

	convert_egv_to_p_block(egv, p, sp, block, nInd, indArr, nIndBath, indBathArr);

	int i, status, **par;
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;
	gsl_vector *x;
	gsl_multimin_function_fdf my_func;

	int *offset, *number, bound[2], resi=NpBlock%size,
	    NpPerSize = NpBlock/size ;
	offset = mkvectori( size+1 );
	number = mkvectori( size );
	if( size > NpBlock ){
		for(i=0; i<NpBlock; i++)
			offset[i] = i;
		for(i=NpBlock; i<size+1; i++)
			offset[i] = NpBlock;
		bound[0] = myrank < NpBlock ? myrank : NpBlock;
		bound[1] = myrank < NpBlock ? myrank+1 : NpBlock;
		//printf("at rank%02d myblock = %d-%d=%d\n", myrank, bound[1], bound[0], bound[1]-bound[0]);
	}
	else{
		for(i=0; i<size; i++)
			offset[i] = i * (NpPerSize);
		for(i=0; i<resi+1; i++)
			offset[i] += i;
		for(i=resi+1; i<size; i++)
			offset[i] += resi;
		offset[i] = NpBlock;
		bound[0] = offset[myrank];
		if( myrank+1 == size ){
			bound[1] = NpBlock;
		}
		else	bound[1] = offset[myrank+1];
	}
	for(i=0; i<size; i++)
		number[i] = offset[i+1] - offset[i];
	//for(i=0; i<size; i++)	printf("rank%d: number[%d] = %d\n", myrank, i, number[i]);

	int npar = 11;
	par = (int **) malloc( npar*sizeof(int*) );
	par[0] = offset;
	par[1] = number;
	par[2] = bound;
	par[3] = &n;

	par[4] = &sp;
	par[5] = &block;
	par[6] = &nInd;
	par[7] = indArr;
	par[8] = &nIndBath;
	par[9] = indBathArr;
	par[10]= &NpBlock;

	my_func.n	= NpBlock;
	my_func.f	= my_f_block_gslpar;
	my_func.df	= my_df_block;
	my_func.fdf	= my_fdf_block;
	my_func.params	= par;

	x = gsl_vector_alloc(NpBlock);
	for(i=0; i<NpBlock; i++)	gsl_vector_set(x, i, p[i]);

	T = gsl_multimin_fdfminimizer_conjugate_fr;
	s = gsl_multimin_fdfminimizer_alloc(T, NpBlock);

	int ex = ( count > 7 ) ? 10 : count+2; 
	double step = pow(10,-ex);
	gsl_multimin_fdfminimizer_set(s, &my_func, x, step, 2e-4);
	*iter = 0;
	int find = 1, allfind;

	do{
		*iter = *iter + 1;
		status = gsl_multimin_fdfminimizer_iterate (s);
		status = gsl_multimin_test_gradient (s->gradient, Converge);
		if( status == GSL_CONTINUE )	find = 0;

		for(i=0; i<NpBlock; i++)	p[i] = gsl_vector_get(s->x, i);
		MPI_Bcast(p, NpBlock, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		for(i=0; i<NpBlock; i++)	gsl_vector_set(s->x, i, p[i]);
		MPI_Allreduce(&find, &allfind, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		if( allfind == size )	break;
		else if( allfind > 0 ){
			if(myrank<1) {
				printf("iter%d rank%d find%d allfind%d\n", *iter, myrank, find, allfind);
				printf("Inconsistency in minimization!\n");	fflush(stdout);
			}
			exit(1);
		}
		else find = 1;
	}while(status == GSL_CONTINUE && *iter < Maxmin);

	for(i=0; i<NpBlock; i++)	p[i] = gsl_vector_get(s->x, i);
	*minimal = s->f;

	gsl_multimin_fdfminimizer_free(s);
	gsl_vector_free(x);

	convert_p_to_egv_block(egv, p, sp, block, nInd, indArr, nIndBath, indBathArr);
	//sort_bath_zeroSOC(&egv, sp);
	char paratemp[1024];
	if( myrank == 0 ){
		sprintf(paratemp, "rank%d minimize %d block %d sp %d (iter %d)", myrank, n, block, sp, *iter);
		print_egv(paratemp, egv);
	}

	free(par);
	free(offset);
	free(number);

}//end of minimize_zeroSOC

void minimize_blockSym(int count, PNSpara egv, double *minimal, int *iter, double Converge, int n, int sp, int block){
	infoHybblock iHyb( "inputs/HYBBLOCK" );
	int NpBlock	= iHyb.getNp( block ) ;			if(myrank<1) cout<< "NpBlock "<<block<<" (symmetrized) : "<< NpBlock << endl ;
	int nInd	= iHyb.blockIndArr[block].size();
	int indArr[tNC];
	for(int i=0; i<nInd; i++ ) indArr[i] = iHyb.blockIndArr[block][i] + sp ;

	int nIndBath	= iHyb.nbathBlockArr[block];
	int indBathArr[Np];
	for(int i=0; i<nIndBath; i++ ) indBathArr[i] = iHyb.bathBlockIndArr[block][i]*(iHyb.ispin+1) + sp ;
	if(myrank<1){
		printf("\tindices of bath-orbital : ");
	       	for(int i=0; i<nIndBath; i++) printf("%d ",indBathArr[i]);
		printf("\n") ;
	}

	//int nIndBathSym = nIndBath/nInd;

	if( NpBlockCap < NpBlock ){ if(myrank<1) printf("Error :: NpBlockCap<NpBlock. %d<%d\n",NpBlockCap,NpBlock); exit(1); }
	double p[NpBlockCap];

	convert_egv_to_p_blockSym(egv, p, sp, block, nInd, indArr, nIndBath, indBathArr);

	int i, status, **par;
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;
	gsl_vector *x;
	gsl_multimin_function_fdf my_func;

	int *offset, *number, bound[2], resi=NpBlock%size,
	    NpPerSize = NpBlock/size ;
	offset = mkvectori( size+1 );
	number = mkvectori( size );
	if( size > NpBlock ){
		for(i=0; i<NpBlock; i++)
			offset[i] = i;
		for(i=NpBlock; i<size+1; i++)
			offset[i] = NpBlock;
		bound[0] = myrank < NpBlock ? myrank : NpBlock;
		bound[1] = myrank < NpBlock ? myrank+1 : NpBlock;
		//printf("at rank%02d myblock = %d-%d=%d\n", myrank, bound[1], bound[0], bound[1]-bound[0]);
	}
	else{
		for(i=0; i<size; i++)
			offset[i] = i * (NpPerSize);
		for(i=0; i<resi+1; i++)
			offset[i] += i;
		for(i=resi+1; i<size; i++)
			offset[i] += resi;
		offset[i] = NpBlock;
		bound[0] = offset[myrank];
		if( myrank+1 == size ){
			bound[1] = NpBlock;
		}
		else	bound[1] = offset[myrank+1];
	}
	for(i=0; i<size; i++)
		number[i] = offset[i+1] - offset[i];
	//for(i=0; i<size; i++)	printf("rank%d: number[%d] = %d\n", myrank, i, number[i]);

	int npar = 11;
	par = (int **) malloc( npar*sizeof(int*) );
	par[0] = offset;
	par[1] = number;
	par[2] = bound;
	par[3] = &n;

	par[4] = &sp;
	par[5] = &block;
	par[6] = &nInd;
	par[7] = indArr;
	par[8] = &nIndBath;
	par[9] = indBathArr;
	par[10]= &NpBlock;

	my_func.n	= NpBlock;
	my_func.f	= my_f_block_gslparSym;
	my_func.df	= my_df_blockSym;
	my_func.fdf	= my_fdf_blockSym;
	my_func.params	= par;

	x = gsl_vector_alloc(NpBlock);
	for(i=0; i<NpBlock; i++)	gsl_vector_set(x, i, p[i]);

	T = gsl_multimin_fdfminimizer_conjugate_fr;
	s = gsl_multimin_fdfminimizer_alloc(T, NpBlock);

	int ex = ( count > 7 ) ? 10 : count+2; 
	double step = pow(10,-ex);
	gsl_multimin_fdfminimizer_set(s, &my_func, x, step, 2e-4);
	*iter = 0;
	int find = 1, allfind;

	do{
		*iter = *iter + 1;
		status = gsl_multimin_fdfminimizer_iterate (s);
		status = gsl_multimin_test_gradient (s->gradient, Converge);
		if( status == GSL_CONTINUE )	find = 0;

		for(i=0; i<NpBlock; i++)	p[i] = gsl_vector_get(s->x, i);
		MPI_Bcast(p, NpBlock, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		for(i=0; i<NpBlock; i++)	gsl_vector_set(s->x, i, p[i]);
		MPI_Allreduce(&find, &allfind, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		if( allfind == size )	break;
		else if( allfind > 0 ){
			if(myrank<1) {
				printf("iter%d rank%d find%d allfind%d\n", *iter, myrank, find, allfind);
				printf("Inconsistency in minimization!\n");	fflush(stdout);
			}
			exit(1);
		}
		else find = 1;
	}while(status == GSL_CONTINUE && *iter < Maxmin);

	for(i=0; i<NpBlock; i++)	p[i] = gsl_vector_get(s->x, i);
	*minimal = s->f;

	gsl_multimin_fdfminimizer_free(s);
	gsl_vector_free(x);

	convert_p_to_egv_blockSym(egv, p, sp, block, nInd, indArr, nIndBath, indBathArr);
	//sort_bath_zeroSOC(&egv, sp);
	char paratemp[1024];
	if( myrank == 0 ){
		sprintf(paratemp, "rank%d minimize %d block %d sp %d (iter %d)", myrank, n, block, sp, *iter);
		print_egv(paratemp, egv);
	}

	free(par);
	free(offset);
	free(number);

}//end of minimize_zeroSOC

void minimize_block(int count, PNSpara egv, double *minimal, int *iter, double Converge, int n, int sp, int block){
	infoHybblock iHyb( "inputs/HYBBLOCK" );
	if( iHyb.blockSymIndArr[block]<1)
		minimize_blockAll(count, egv, minimal, iter, Converge, n, sp, block);
	else
		minimize_blockSym(count, egv, minimal, iter, Converge, n, sp, block);
}
