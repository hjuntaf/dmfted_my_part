#include <cstdlib>
#include <cmath>
#include <string>
#define  SQRTTWO 1.414213562373095
#define hSQRTTWO 0.707106781186547

class comp {
	public :
		int i,j ;
		double dat ;
		comp(void) {}
		comp( int ii, int jj, double dd ) {
			i = ii ; j = jj ; dat = dd ;
		}
		void setDat( int ii, int jj, double dd ) {
			i = ii ; j = jj ; dat = dd ;
		}
} ;
class compZ {
	public :
		int i,j ;
		gsl_complex dat ;
		compZ(void) {}
		compZ( int ii, int jj, double ddre, double ddim ) {
			i = ii ; j = jj ; dat = gsl_complex_rect(ddre,ddim) ;
		}
		void setDat( int ii, int jj, double ddre, double ddim ) {
			i = ii ; j = jj ; dat = gsl_complex_rect(ddre,ddim) ;
		}
} ;
class Op {
	public :
		int dim , ndat ; 
		comp *cDat ; 
		Op(void) {}
		Op( char opname[] ) {
			std::string opnstr (opname) ;
			if( opnstr.find("TrevUnit") < std::string::npos ){
				initTrev(0) ;
			}
			else if( opnstr.find("Trev") < std::string::npos ){
				initTrev(1) ;
			}
			else {
				printf("ERROR :: Constructor() of Op : invalid 'opname[]' .\n" ) ;
			}
		}
		Op( int dd , int nn ) {
			dim  = dd ;
			ndat = nn ;
			cDat = (comp*)malloc( sizeof(comp)*ndat ) ;
		}
		void mkOp( int dd , int nn ) {
			dim  = dd ;
			ndat = nn ;
			cDat = (comp*)malloc( sizeof(comp)*ndat ) ;
		}
		void freeOp() {
			free(cDat) ;
		}
		void showOp() {
			for( int a ; a<ndat ; a++ ) {
				printf("%4d %4d %19.16f\n", cDat[a].i, cDat[a].j, cDat[a].dat ) ;
			}
		}
		void initTrev(int fermionsgn) {
			mkOp(6,6);
			if(fermionsgn) {
			cDat[0].setDat( 3, 0,  1 );
			cDat[1].setDat( 2, 1, -1 );
			cDat[2].setDat( 1, 2,  1 );
			cDat[3].setDat( 0, 3, -1 );
			cDat[4].setDat( 5, 4,  1 );
			cDat[5].setDat( 4, 5, -1 );
			}
			else {
			cDat[0].setDat( 3, 0,  1 );
			cDat[1].setDat( 2, 1,  1 );
			cDat[2].setDat( 1, 2,  1 );
			cDat[3].setDat( 0, 3,  1 );
			cDat[4].setDat( 5, 4,  1 );
			cDat[5].setDat( 4, 5,  1 );
			}
		}
} ;
class OpZ {
	public :
		int dim , ndat ; 
		compZ *cDat ; 
		OpZ(void) {}
		OpZ( int dd , int nn ) {
			dim  = dd ;
			ndat = nn ;
			cDat = (compZ*)malloc( sizeof(compZ)*ndat ) ;
		}
		void mkOp( int dd , int nn ) {
			dim  = dd ;
			ndat = nn ;
			cDat = (compZ*)malloc( sizeof(compZ)*ndat ) ;
		}
		void freeOp() {
			free(cDat) ;
		}
		void showOp() {
			for( int a ; a<ndat ; a++ ) {
				printf("%4d %4d %19.16f %19.16f\n", cDat[a].i, cDat[a].j, GSL_REAL(cDat[a].dat), GSL_IMAG(cDat[a].dat) ) ;
			}
		}
} ;


class Vecop {
	public : 
		Op zj, pj, mj ;
		Vecop(void) {}
		Vecop( char opname[] ) {
			std::string opnstr (opname) ;
			if( opnstr.find("Jop") < std::string::npos ){
				initJ() ;
			}
			else if( opnstr.find("Jaop") < std::string::npos ){
				initJa() ;
			}
			else if( opnstr.find("Lop") < std::string::npos ){
				initL() ;
			}
			else if( opnstr.find("Sop") < std::string::npos ){
				initS() ;
			}
			else if( opnstr.find("Mop") < std::string::npos ){
				initM() ;
			}
			else {
				printf("ERROR :: Constructor() of Vecop : invalid 'opname[]' .\n" ) ;
			}
		}
		void initJ(void) {
			zj.mkOp(6,6) ;
			pj.mkOp(6,4) ;
			mj.mkOp(6,4) ; 
			zj.cDat[0].setDat( 0, 0,  1.5 ) ;
			zj.cDat[1].setDat( 1, 1,  0.5 ) ;
			zj.cDat[2].setDat( 2, 2, -0.5 ) ;
			zj.cDat[3].setDat( 3, 3, -1.5 ) ;
			zj.cDat[4].setDat( 4, 4,  0.5 ) ;
			zj.cDat[5].setDat( 5, 5, -0.5 ) ;

			pj.cDat[0].setDat( 0, 1,  sqrt(3)	) ;
			pj.cDat[1].setDat( 1, 2,  2      	) ;
			pj.cDat[2].setDat( 2, 3,  sqrt(3)	) ;
			pj.cDat[3].setDat( 4, 5,  1      	) ;

			mj.cDat[0].setDat( 1, 0,  sqrt(3)	) ;
			mj.cDat[1].setDat( 2, 1,  2      	) ;
			mj.cDat[2].setDat( 3, 2,  sqrt(3)	) ;
			mj.cDat[3].setDat( 5, 4,  1      	) ;
		}
		void initJa(void) {
			zj.mkOp(6,10) ;
			pj.mkOp(6,8) ;
			mj.mkOp(6,8) ; 
			zj.cDat[0].setDat( 0, 0, -0.5		) ;
			zj.cDat[1].setDat( 1, 1, -1./6		) ;
			zj.cDat[2].setDat( 2, 2,  1./6		) ;
			zj.cDat[3].setDat( 3, 3,  0.5		) ;
			zj.cDat[4].setDat( 4, 4, -5./6		) ;
			zj.cDat[5].setDat( 5, 5,  5./6		) ;
			double dum = 2.*sqrt(2)/3 ;
			zj.cDat[6].setDat( 1, 4, -dum		) ;
			zj.cDat[7].setDat( 4, 1, -dum		) ;
			zj.cDat[8].setDat( 2, 5, -dum		) ;
			zj.cDat[9].setDat( 5, 2, -dum		) ;

			pj.cDat[0].setDat( 0, 1,  1./sqrt(3)	) ;
			pj.cDat[1].setDat( 1, 2, -2./3		) ;
			pj.cDat[2].setDat( 2, 3, -1./sqrt(3)	) ;
			pj.cDat[3].setDat( 4, 5, -5./3		) ;
			pj.cDat[4].setDat( 0, 4,  dum*sqrt(3)	) ;
			pj.cDat[5].setDat( 1, 5,  dum		) ;
			pj.cDat[6].setDat( 4, 2, -dum		) ;
			pj.cDat[7].setDat( 5, 3, -dum*sqrt(3)	) ;

			mj.cDat[0].setDat( 1, 0, -1./sqrt(3)	) ;
			mj.cDat[1].setDat( 2, 1, -2./3		) ;
			mj.cDat[2].setDat( 3, 2, -1./sqrt(3)	) ;
			mj.cDat[3].setDat( 5, 4, -5./3		) ;
			mj.cDat[4].setDat( 4, 0,  dum*sqrt(3)	) ;
			mj.cDat[5].setDat( 5, 1,  dum		) ;
			mj.cDat[6].setDat( 2, 4, -dum		) ;
			mj.cDat[7].setDat( 3, 5, -dum*sqrt(3)	) ;
		}
		void initL(void) {
			zj.mkOp(6,10) ;
			pj.mkOp(6,8) ;
			mj.mkOp(6,8) ; 
			zj.cDat[0].setDat( 0, 0,  1.		) ;
			zj.cDat[1].setDat( 1, 1,  1./3		) ;
			zj.cDat[2].setDat( 2, 2, -1./3		) ;
			zj.cDat[3].setDat( 3, 3, -1.		) ;
			zj.cDat[4].setDat( 4, 4,  2./3		) ;
			zj.cDat[5].setDat( 5, 5, -2./3		) ;
			double dum = sqrt(2)/3 ;
			zj.cDat[6].setDat( 1, 4,  dum		) ;
			zj.cDat[7].setDat( 4, 1,  dum		) ;
			zj.cDat[8].setDat( 2, 5,  dum		) ;
			zj.cDat[9].setDat( 5, 2,  dum		) ;

			pj.cDat[0].setDat( 0, 1,  2./sqrt(3)	) ;
			pj.cDat[1].setDat( 1, 2,  4./3		) ;
			pj.cDat[2].setDat( 2, 3,  2./sqrt(3)	) ;
			pj.cDat[3].setDat( 4, 5,  4./3		) ;
			pj.cDat[4].setDat( 0, 4, -dum*sqrt(3)	) ;
			pj.cDat[5].setDat( 1, 5, -dum		) ;
			pj.cDat[6].setDat( 4, 2,  dum		) ;
			pj.cDat[7].setDat( 5, 3,  dum*sqrt(3)	) ;

			mj.cDat[0].setDat( 1, 0,  2./sqrt(3)	) ;
			mj.cDat[1].setDat( 2, 1,  4./3		) ;
			mj.cDat[2].setDat( 3, 2,  2./sqrt(3)	) ;
			mj.cDat[3].setDat( 5, 4,  4./3		) ;
			mj.cDat[4].setDat( 4, 0, -dum*sqrt(3)	) ;
			mj.cDat[5].setDat( 5, 1, -dum		) ;
			mj.cDat[6].setDat( 2, 4,  dum		) ;
			mj.cDat[7].setDat( 3, 5,  dum*sqrt(3)	) ;
		}
		void initS(void) {
			zj.mkOp(6,10) ;
			pj.mkOp(6,8) ;
			mj.mkOp(6,8) ; 
			zj.cDat[0].setDat( 0, 0,  0.5		) ;
			zj.cDat[1].setDat( 1, 1,  1./6		) ;
			zj.cDat[2].setDat( 2, 2, -1./6		) ;
			zj.cDat[3].setDat( 3, 3, -0.5		) ;
			zj.cDat[4].setDat( 4, 4, -1./6		) ;
			zj.cDat[5].setDat( 5, 5,  1./6		) ;
			double dum = sqrt(2)/3 ;
			zj.cDat[6].setDat( 1, 4, -dum		) ;
			zj.cDat[7].setDat( 4, 1, -dum		) ;
			zj.cDat[8].setDat( 2, 5, -dum		) ;
			zj.cDat[9].setDat( 5, 2, -dum		) ;

			pj.cDat[0].setDat( 0, 1,  1./sqrt(3)	) ;
			pj.cDat[1].setDat( 1, 2,  2./3		) ;
			pj.cDat[2].setDat( 2, 3,  1./sqrt(3)	) ;
			pj.cDat[3].setDat( 4, 5, -1./3		) ;
			pj.cDat[4].setDat( 0, 4,  dum*sqrt(3)	) ;
			pj.cDat[5].setDat( 1, 5,  dum		) ;
			pj.cDat[6].setDat( 4, 2, -dum		) ;
			pj.cDat[7].setDat( 5, 3, -dum*sqrt(3)	) ;

			mj.cDat[0].setDat( 1, 0,  1./sqrt(3)	) ;
			mj.cDat[1].setDat( 2, 1,  2./3		) ;
			mj.cDat[2].setDat( 3, 2,  1./sqrt(3)	) ;
			mj.cDat[3].setDat( 5, 4, -1./3		) ;
			mj.cDat[4].setDat( 4, 0,  dum*sqrt(3)	) ;
			mj.cDat[5].setDat( 5, 1,  dum		) ;
			mj.cDat[6].setDat( 2, 4, -dum		) ;
			mj.cDat[7].setDat( 3, 5, -dum*sqrt(3)	) ;
		}
		void initM(void) {
			zj.mkOp(6,6) ;
			pj.mkOp(6,5) ;
			mj.mkOp(6,5) ; 
			double dum = sqrt(2) ;
			zj.cDat[0].setDat( 4, 4, -1.		) ;
			zj.cDat[1].setDat( 5, 5,  1.		) ;
			zj.cDat[2].setDat( 1, 4, -dum		) ;
			zj.cDat[3].setDat( 4, 1, -dum		) ;
			zj.cDat[4].setDat( 2, 5, -dum		) ;
			zj.cDat[5].setDat( 5, 2, -dum		) ;

			pj.cDat[0].setDat( 4, 5, -2.		) ;
			pj.cDat[1].setDat( 0, 4,  dum*sqrt(3)	) ;
			pj.cDat[2].setDat( 1, 5,  dum		) ;
			pj.cDat[3].setDat( 4, 2, -dum		) ;
			pj.cDat[4].setDat( 5, 3, -dum*sqrt(3)	) ;

			mj.cDat[0].setDat( 5, 4, -2.		) ;
			mj.cDat[1].setDat( 4, 0,  dum*sqrt(3)	) ;
			mj.cDat[2].setDat( 5, 1,  dum		) ;
			mj.cDat[3].setDat( 2, 4, -dum		) ;
			mj.cDat[4].setDat( 3, 5, -dum*sqrt(3)	) ;
		}
		void freeVecop(void) {
			zj.freeOp() ;
			pj.freeOp() ;
			mj.freeOp() ;
		}
		void showVop(void) {
			printf("zj :\n") ; zj.showOp() ;
			printf("pj :\n") ; pj.showOp() ;
			printf("mj :\n") ; mj.showOp() ;
		}
} ;

class VecopZ {
	public : 
		OpZ *op ;
		Op zj, pj, mj ;
		int nop ;
		VecopZ(void) {
			nop=6;
			init() ;
		}
		VecopZ( char opname[] ) {
			std::string opnstr (opname) ;
			if( opnstr.find("D4h") < std::string::npos ){
				nop=5;
				init() ;
			}
			if( opnstr.find("Tt2gjeff") < std::string::npos ){
				nop=1;
				initTt2gjeff() ;
			}
			else {
				printf("ERROR :: Constructor() of VecopZ : invalid 'opname[]' .\n" ) ;
			}
		}
		void init(void) {
			op = (OpZ *)malloc( sizeof(OpZ)*nop );

			//Time-reversal
			op[0].mkOp(6,6) ;
			op[0].cDat[0].setDat( 0, 3, -1,0 );
			op[0].cDat[1].setDat( 1, 2,  1,0 );
			op[0].cDat[2].setDat( 2, 1, -1,0 );
			op[0].cDat[3].setDat( 3, 0,  1,0 );
			op[0].cDat[4].setDat( 4, 5, -1,0 );
			op[0].cDat[5].setDat( 5, 4,  1,0 );

			//C4
			op[1].mkOp(6,6) ;
			op[1].cDat[0].setDat( 0, 0,  hSQRTTWO, hSQRTTWO );
			op[1].cDat[1].setDat( 1, 1, -hSQRTTWO, hSQRTTWO );
			op[1].cDat[2].setDat( 2, 2, -hSQRTTWO,-hSQRTTWO );
			op[1].cDat[3].setDat( 3, 3,  hSQRTTWO,-hSQRTTWO );
			op[1].cDat[4].setDat( 4, 4, -hSQRTTWO, hSQRTTWO );
			op[1].cDat[5].setDat( 5, 5, -hSQRTTWO,-hSQRTTWO );

			//C2
			op[2].mkOp(6,6) ;
			op[2].cDat[0].setDat( 0, 0,  0., 1 );
			op[2].cDat[1].setDat( 1, 1,  0.,-1 );
			op[2].cDat[2].setDat( 2, 2,  0., 1 );
			op[2].cDat[3].setDat( 3, 3,  0.,-1 );
			op[2].cDat[4].setDat( 4, 4,  0.,-1 );
			op[2].cDat[5].setDat( 5, 5,  0., 1 );

			//C2p
			op[3].mkOp(6,10) ;
			op[3].cDat[0].setDat( 0, 3,  hSQRTTWO,-hSQRTTWO );
			op[3].cDat[1].setDat( 1, 2, -hSQRTTWO, hSQRTTWO/3. );
			op[3].cDat[2].setDat( 2, 1,  hSQRTTWO, hSQRTTWO/3. );
			op[3].cDat[3].setDat( 3, 0, -hSQRTTWO,-hSQRTTWO );

			op[3].cDat[4].setDat( 1, 5,        0., 2./3. );
			op[3].cDat[5].setDat( 2, 4,        0.,-2./3. );
			op[3].cDat[6].setDat( 4, 2,        0.,-2./3. );
			op[3].cDat[7].setDat( 5, 1,        0., 2./3. );

			op[3].cDat[8].setDat( 4, 5,  hSQRTTWO, SQRTTWO/3. );
			op[3].cDat[9].setDat( 5, 4, -hSQRTTWO, SQRTTWO/3. );

			//C2pp
			op[4].mkOp(6,10) ;
			op[4].cDat[0].setDat( 0, 3,        -1,        0 );
			op[4].cDat[1].setDat( 1, 2,    -1./3.,    -2./3 );
			op[4].cDat[2].setDat( 2, 1,     1./3.,    -2./3 );
			op[4].cDat[3].setDat( 3, 0,         1,        0 );

			op[4].cDat[4].setDat( 1, 5,  SQRTTWO/3., -SQRTTWO/3. );
			op[4].cDat[5].setDat( 2, 4,  SQRTTWO/3.,  SQRTTWO/3. );
			op[4].cDat[6].setDat( 4, 2, -SQRTTWO/3.,  SQRTTWO/3. );
			op[4].cDat[7].setDat( 5, 1, -SQRTTWO/3., -SQRTTWO/3. );

			op[4].cDat[8].setDat( 4, 5,      2./3,     1./3 );
			op[4].cDat[9].setDat( 5, 4,     -2./3,     1./3 );

			//R(2*pi)
			op[5].mkOp(6,6) ;
			op[5].cDat[0].setDat( 0, 0,  -1., 0 );
			op[5].cDat[1].setDat( 1, 1,  -1., 0 );
			op[5].cDat[2].setDat( 2, 2,  -1., 0 );
			op[5].cDat[3].setDat( 3, 3,  -1., 0 );
			op[5].cDat[4].setDat( 4, 4,  -1., 0 );
			op[5].cDat[5].setDat( 5, 5,  -1., 0 );
		}
		void initTt2gjeff(void) {
			op = (OpZ *)malloc( sizeof(OpZ)*nop );
			int jj=0;

			//C4
			op[0].mkOp(6,16) ;
			op[0].cDat[jj].setDat( 0, 2,-1./SQRTTWO , 0.         ); jj++;
			op[0].cDat[jj].setDat( 0, 4, 0.         ,-1./SQRTTWO ); jj++;
			op[0].cDat[jj].setDat( 1, 0, sqrt(2./3.), 0.         ); jj++;
			op[0].cDat[jj].setDat( 1, 3,-sqrt(1./6.), 0.         ); jj++;
			op[0].cDat[jj].setDat( 1, 5, 0.         ,-sqrt(1./6.)); jj++;
			op[0].cDat[jj].setDat( 2, 1, sqrt(2./3.), 0.         ); jj++;
			op[0].cDat[jj].setDat( 2, 2, sqrt(1./6.), 0.         ); jj++;
			op[0].cDat[jj].setDat( 2, 4, 0.         ,-sqrt(1./6.)); jj++;
			op[0].cDat[jj].setDat( 3, 2, 1./SQRTTWO , 0.         ); jj++;
			op[0].cDat[jj].setDat( 3, 4, 0.         ,-1./SQRTTWO ); jj++;
			op[0].cDat[jj].setDat( 4, 0,-sqrt(1./3.), 0.         ); jj++;
			op[0].cDat[jj].setDat( 4, 3,-sqrt(1./3.), 0.         ); jj++;
			op[0].cDat[jj].setDat( 4, 5, 0.         ,-sqrt(1./3.)); jj++;
			op[0].cDat[jj].setDat( 5, 1, sqrt(1./3.), 0.         ); jj++;
			op[0].cDat[jj].setDat( 5, 2,-sqrt(1./3.), 0.         ); jj++;
			op[0].cDat[jj].setDat( 5, 4, 0.         , sqrt(1./3.)); jj++;
		}
		void freeVecop(void) {
			for( int a=0 ; a<nop ; a++ ) 
				op[a].freeOp() ;
		}
		void showVop(void) {
			for( int a=0 ; a<nop ; a++ ) {
				printf("op%d :\n",a) ; op[a].showOp() ;
			}
		}
} ;
