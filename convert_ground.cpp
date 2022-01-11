int search_gindex(sectorinfo *ginfo, int refer, int block, int Nblocks){
	int i;
	for(i=0; i<Nblocks; i++){
		if( ginfo[i].refer == refer && ginfo[i].block == block )	return i;
	}
	my_error("Error in search_gindex!");

	return 0;
}

gsl_complex* convert_ground(double *ginformation, gsl_complex *ground, int *gind, int *gnd_conv, int *gndblock_conv ){
	int Blocks_origin = (Ns+1)*(Ns+1);
	typebasis **basis_orig, **basis_converted;
    int Powns2dum    = (int)pow(2,2*Ns);
	basis_orig       = mkmatrixb(Powns2dum, 2);
	basis_converted  = mkmatrixb(Powns2dum, 2);
	sectorinfo *ginfo_orig, *ginfo_converted;

	ginfo_orig = (sectorinfo *) malloc( Blocks_origin*sizeof(sectorinfo) );
	build_bitbasis(Ns, basis_orig, ginfo_orig, 0);

	ginfo_converted = (sectorinfo *) malloc( Blocks*sizeof(sectorinfo) );
	build_bitbasis(Ns, basis_converted, ginfo_converted, 1);

	int gnd, gndblock, gindex, gnd_converted, gndblock_converted;
	gnd = (int) ginformation[1];
	gndblock = (int) ginformation[2];

	gindex = search_gindex(ginfo_orig, gnd, gndblock, Blocks_origin); 
	int new_gindex = ginfo_orig[gindex].ptl;

	//printf("Original expect_hop :: \n");
	//for( int mu=0 ; mu<tNC ; mu++ ) {
	//	double dum = expect_hop( basis, &ginfo[gindex], ground, mu, mu );
	//	printf("%d,%d : %19.16f\n", mu, mu, dum );
	//}
	//for( int mu=0 ; mu<NC ; mu++ ) for( int nu=0 ; nu<NC ; nu++ ) {
	//	double dum = expect_hop( basis, &ginfo[gindex], ground, 2*mu, 2*nu );
	//	printf("%d,%d : %19.16f\n", 2*mu, 2*nu, dum );
	//	dum = expect_hop( basis, &ginfo[gindex], ground, 2*mu+1, 2*nu+1 );
	//	printf("%d,%d : %19.16f\n", 2*mu+1, 2*nu+1, dum );
	//}

	printf("convert_ground:: gindex = %d, new_gindex = %d\n", gindex, new_gindex);
	gnd_converted		= ginfo_converted[new_gindex].refer;
	gndblock_converted	= ginfo_converted[new_gindex].block;
	if( ginfo_orig[gindex].ptl != ginfo_converted[new_gindex].ptl )	my_error("wrong conversion error in convert_ground");

	gsl_complex *ground_new = mkgscvectord( gndblock_converted );
	int *table = mkvectori( Powns2dum );
#pragma omp parallel for default(shared)
	for(int i=0; i<gndblock_converted; i++)	ground_new[i] = zero;
#pragma omp parallel for default(shared)
	for(int i=0; i<Powns2dum; i++) table[i] = 0;
#pragma omp parallel for default(shared)
	for(int i=0; i<Powns2dum; i++){
		int index = (basis_converted[i][0]<<Ns) + basis_converted[i][1];
		table[index] = i;
	}
	typebasis bin[2];
	for(int i=gnd; i<gnd+gndblock; i++){
		bin[0] = basis_orig[i][0];
		bin[1] = basis_orig[i][1];
		int index = table[(bin[0]<<Ns) + bin[1]];
		GSL_REAL(ground_new[index-gnd_converted]) = GSL_REAL(ground[i-gnd]);
		GSL_IMAG(ground_new[index-gnd_converted]) = GSL_IMAG(ground[i-gnd]);
	}
	*gind = new_gindex;
	*gnd_conv = gnd_converted;
	*gndblock_conv = gndblock_converted;
	free(ground);

	printf("Transforming to Tt2gjeff\n"); fflush(stdout);
	gsl_complex *ground_newT= mkgscvectord( gndblock_converted );
	obtainGroundT( ground_new, ground_newT, gnd_converted, gndblock_converted, basis_converted, table);
	free(ground_new);

	free(table);
	free(ginfo_orig);
	free(ginfo_converted);
	freematrixb(basis_orig, Powns2dum);
	freematrixb(basis_converted, Powns2dum);
	return ground_newT;
}

