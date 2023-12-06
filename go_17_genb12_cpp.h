
void GEN_BANDES_12::GetStartB2(int ip) {//set  rows 3_9 column 1
	char const *tp[20] = {// 3 out of 6 ordered 
		"012345", "345012", "013245", "245013", "014235", "235014", "015234", "234015",
		"023145", "145023", "024135", "135024", "025134", "134025",
		"034125", "125034", "035124", "124035", "045123", "123045"
	};
		char tpw[7];	strcpy( tpw , tp[ip]);
	for (int i = 0; i < 6; i++)boxd[i] = 0x1ff;
	for (int j = 0, jc = 27; j < 6; j++, jc += 9) {
		int ic = tpw[j] - '0', c0 = tc[ic], bit = 1 << c0;
		grid0[jc] = c0;		zsol[jc] = c0+'1';
		rowd[j] = 0x1ff ^ bit;
		if (j < 3)boxd[0] ^= bit; else  boxd[3] ^= bit;
	}
}
void GEN_BANDES_12::Start(int mode) {
	modeb12 = mode;
	myband1.Initstd();
	zsol[81] = 0;	nb12 = 0;
}
void GEN_BANDES_12::NewBand1(int iw) {
	go_back = 0;
	i1t16 = iw;	it16 = tn6_to_416[iw];
	ib1check = i1t16;
	myband1.InitG12(it16);
	memcpy(grid0, myband1.band0, sizeof myband1.band0);
	memcpy(gcheck, myband1.band0, sizeof myband1.band0);
	n_auto_b1 = tblnauto[it16]; //ia 0-415 not index
	t_auto_b1 = &automorphsp[tblnautostart[it16]];

	strcpy(zsol, myband1.band);
	for (int i = 0; i < 9; i++) // init columns status
		cold[i] = 0x1ff ^ myband1.gangster[i];
	memcpy(coldf, cold, sizeof coldf);
	zsol[27] = 0;
	if(op.ton)	cout << myband1.band << "i1t16=" << i1t16 << " it16=" 
		<< it16		<< " n auto morphs=" << n_auto_b1 << endl;

	ntc = 0;
	BitsInTable32(tc, ntc, cold[0]);// first col 6 digits in table
	for (ip20 = 0; ip20 < 20; ip20++) {//0;ip<20 setup initial values for rows columns
		GetStartB2(ip20);
		Find_band2B();
		if (go_back)return;
	}
}
int GEN_BANDES_12::Band2Check() {
	int * zs0 = &gcheck[27];
	n_auto_b1b2 = 0;
	if (n_auto_b1) {
		for (int imorph = 0; imorph < n_auto_b1; imorph++) {
			BANDMINLEX::PERM &p = t_auto_b1[imorph];
			int band[27];// morph the band
			for (int i = 0; i < 9; i++) {
				band[i] = p.map[zs0[p.cols[i]]];
				band[i + 9] = p.map[zs0[p.cols[i] + 9]];
				band[i + 18] = p.map[zs0[p.cols[i] + 18]];
			}
			int ir = G17ComparedOrderedBand(zs0, band);
			if (ir == 1)				return 1;
			else if (!ir) {// auto morph b1 b2 store it for later
				t_auto_b1b2[n_auto_b1b2++] = p;
			}
		}
	}
	n_auto_b2b1 = 0;// possible automorph after perm b1b2
	if (i1t16 == ib2check) {// must try perm bands 12 auto morphs
		int b23[3][9];
		for (int i = 0; i < 3; i++) {// morph band1 to band2 minlex
			register int * rrd = b23[i], *rro = &gcheck[9 * i];
			for (int j = 0; j < 9; j++)
				rrd[j] = pcheck2.map[rro[pcheck2.cols[j]]];
		}
		int ir = G17ComparedOrderedBand(zs0, b23[0]);// is it same as base
		if (ir == 1) 			return 1;
		else if (!ir)// auto morph b1 b2 store it for later
			t_auto_b2b1[n_auto_b2b1++].InitBase(ib2check);
		// must also test all auto morphs b2b1
		for (int imorph = 0; imorph < n_auto_b1; imorph++) {// same automorphs b1 b2
			BANDMINLEX::PERM &pp = t_auto_b1[imorph];
			int b23_a[3][9];
			for (int i = 0; i < 3; i++) {
				register int * rrd = b23_a[i], *rro = b23[i];
				for (int j = 0; j < 9; j++)		rrd[j] = pp.map[rro[pp.cols[j]]];
			}
			int ir = G17ComparedOrderedBand(zs0, b23_a[0]);
			if (ir == 1)return 1;
			else if (!ir)// auto morph b1 b2 store it for later
				t_auto_b2b1[n_auto_b2b1++] = pp;
		}
	}
	return 0;
}

int GEN_BANDES_12::F17Novalid1_2() {
	if (!op.t18) {
		int lim = (op.p1) ? 5 : 6;
		if (t16_min_clues[myband1.i416] == 6)
			if (t16_min_clues[myband2.i416] >= lim) {
				p_cpt2g[9] ++;
				//cout << " bands 1+2 with no valid solution "
				//	<< myband1.i416 << " " << myband2.i416 << " " << endl;
				return 1;
			}
	}
	if (op.t18 && op.b3low) return 0;
	//if (op.b2slice) {
		//int ix = t416_to_n6[it16_2];
		//if (ix < op.b2_is) return 1;
		//if (ix > op.b2) return 1;
	//}
	//else 
	if (op.b2<416) {
		if( t416_to_n6[it16_2] != op.b2) return 1;
		if (op.b2start) {
			char* wc = op.b2start;
			int n = (int)strlen(wc);
			if(strncmp(myband2.band,wc,n)) return 1; 
		}
	}
	return 0;
}
int diagbugb2(int * g) {
			  //214865937768932415935714268391548672546271893872693541
	char* zt = "214865937768932415935714268391548672546271893872693541";
	for (int i = 27; i < 54; i++) 	if ((zt[i - 27] - '1') != g[i]) return 0;
	return 1;
}

void GEN_BANDES_12::Find_band2B() {
	int * zs0= &grid0[27];
	register int  *crcb, bit;
	int cd[9], rd[6], bd[6];
	memcpy(rd, rowd, sizeof rd);
	memcpy(cd, coldf, sizeof cd);
	memcpy(bd, boxd, sizeof bd);
	char * zs = zsol;
	// now loop over the 24 cells not yet assigned in the band to fill the band (use relative cell) 
	int ii = -1, free[24];
	uint32_t d;
nextii:
	ii++;
	{	crcb = tgen_band_cat[ii];//cell_r_c_b  24 cells to fill
		register int fr0 = cd[crcb[2]] & bd[crcb[3]], fr = rd[crcb[1]] & fr0;
		if (crcb[4])if (_popcnt32(fr0) < 3) goto back; // 3 clues needed here
		if (!fr)goto back;
		free[ii] = fr;
	}
	goto next_first;
next:// erase previous fill and look for next
	crcb = tgen_band_cat[ii];
	d = zs0[crcb[0]];
	bit = 1 << d;
	rd[crcb[1]] ^= bit; cd[crcb[2]] ^= bit; bd[crcb[3]] ^= bit;
	if (!free[ii])goto back;
	{
	next_first:
		crcb = tgen_band_cat[ii];// be sure to have the good one
		bitscanforward(d, free[ii]);
		bit = 1 << d;
		free[ii] ^= bit;
		zs[crcb[0] + 27] = (char)(d + '1');
		zs0[crcb[0]] = d;
		rd[crcb[1]] ^= bit; cd[crcb[2]] ^= bit; bd[crcb[3]] ^= bit;
		if (ii < 23) goto nextii;
		// this is a valid band, check if lexically minimale 
		int ir = bandminlex.Getmin(zs0, &pband2, 0);
		if (ir < 0) return; //would be bug  did not come in enumeration
		it16_2 = pband2.i416;
		i2t16 = t416_to_n6[it16_2];
		if (i2t16 < i1t16)goto next;// not canonical
		wb1b2.ntab1b2 = wb2b1.ntab1b2 = 0;
		if (n_auto_b1) {
			wb1b2.Init(grid0);
			if (wb1b2.Check(*this) == 1)goto next;
			if (op.p2 && (!op.p2b)) {// compute an apply auto morphs 
				if (i2t16 == i1t16) {
					int ir = wb2b1.Initb2b1(*this);
					if (ir == 1)goto next;
					if (!ir)wb2b1.tab1b2[wb2b1.ntab1b2++].InitBase(i2t16);
					if (wb2b1.Check(*this) == 1)goto next;
				}
			}
		}
		nb12++;
		if (ValidBand2()) { cout << "stop b2" << endl;	go_back = 1; return; }
		goto next;
	}
back:
	if (--ii >= 0) goto next;
}

void GEN_BANDES_12::ValidInitGang() {
	for (int i = 0; i < 9; i++) {// init columns status
		cold[i] = 0x1ff;
		for (int j = 0; j < 6; j++)	cold[i] ^= 1 << grid0[i + 9 * j];
		gangb12[i] = 0x1ff ^ cold[i];
	}
	memcpy(gangcols, cold, sizeof gangcols);
}

void Process81_Split() {
	const char* thexa = { "0123456789ABCDEF" };
	if (genb12.nb12 & 63) return;
	uint64_t sliceclosed = (genb12.nb12 >> 6) - 1;
	GINT64 w64; w64.u64 = 0;
	uint64_t i8deb = p_cptg[1] >> 3;
	uint32_t shift = (uint32_t)(p_cptg[1] - (i8deb << 3));
	int nn = (uint32_t)(p_cpt[1] - p_cptg[1]);
	cout <<sliceclosed<< " slice \t" << p_cptg[1] << " \t" << p_cpt[1]
		<< "\t8d=" << i8deb << "\tshift=" << shift << " nn=" << nn << endl;
	char wss[11];
	sprintf(wss,"%8d ;",(int) sliceclosed);
	wss[10] = 0;
	fout1 << wss;
	p_cptg[1] = p_cpt[1];
	register uint32_t r1 , r2 , rshift = shift;
	register uint8_t* p = bitfield_sgs + i8deb;
	r1 = *p++;
	r1 >>=  shift;// next byte n high bits
	nn -= 8-shift; 
	int nout = 40;
	while (nn>0) {
		r2 = *p++;
		r2 <<=8- shift;
		r1 |= r2;
		register char byte = r1;
		fout1 << thexa[byte & 15] << thexa[(byte >> 4) & 15];
		r1>>= 8;
		nn -= 8;
		nout--;
		if (!nout) {
			fout1<<endl << wss;
			nout = 40;
		}
	}
	register char byte = r1;
	fout1 << thexa[byte & 15] << thexa[(byte >> 4) & 15] << endl;

}
int GEN_BANDES_12::ValidBand2() {
	if (g17b.aigstop)return 1;
	myband2.InitBand2_3(it16_2, &zsol[27], pband2);
	//_______________________ std process
	if (modeb12 < 10) {
		nband3 = 0;
		if((nb12 >> 6) < op.first) return 0;// here restart value, kept untouched if no band 3 found
		{// print a restart point every 64 bands 1+2 seen
			uint64_t w = nb12-1, w1 = w >> 6;
			w &= 63;
			if (w == 0) {
				long tfin = GetTimeMillis();
				if ((w1 > op.last)) return 1;
				if(0)
				cout << "new slice =\t" << w1 << "\tmil=" << (tfin - sgo.tdeb) / 1000
					<< "\tnb2=" << p_cpt2g[0] << "\tnb3=" << p_cpt2g[1] << endl;
			}
		}
		if ((nb12 >> 6) > op.last) return 1;
		ValidInitGang();// also called from existing 17 in test
		if(F17Novalid1_2())return 0;
		Find_band3B();
		return 0;
	}
	//______________________ testing options 
	if (modeb12 ==11) {	// enumeration test
		for (int i = 0; i < 9; i++) {// init columns status
			cold[i] = 0x1ff;
			for (int j = 0; j < 6; j++)	cold[i] ^= 1 << grid0[i + 9 * j];
		}
		memcpy(gangcols, cold, sizeof gangcols);
		if (op.b2<416 && (t416_to_n6[it16_2] != op.b2)) return 0;

		Find_band3B(0);
		if (nband3) {	p_cpt[0]++;	p_cpt[1] += nband3;	}
		if (((nb12 >> 6) > op.last)) return 1;
	}
	return 0;
}
void GEN_BANDES_12::OutEntry() {
	p_cpt2g[0]++;
	p_cpt2g[1] += nband3;   // update bands3 count
	if (op.out_entry < 0) return;
	char ws[21];
	for (int i = 0; i < nband3; i++) {
		fout1 << myband1.band << myband2.band
			<< bands3[i].band;
		if (op.ton) {
			sprintf(ws,";%3d;%3d;%3d;%5d\n", i1t16, i2t16, t416_to_n6[bands3[i].i416], int(nb12 >> 6));
			fout1 << ws;
		}
		else fout1 << endl;
	}
}

void GEN_BANDES_12::Find_band3B(int m10) {
	//if (sgo.bfx[6]) {
		//cout << myband2.band << " new band 2" << endl;
	//}
	register int* crcb, bit;
	nband3 =  nsgchecked=0;
	int* rd = rowdb3, * cd = cold, * bd = boxdb3; // to updates rows cols boxes
	char* zs = zsol;
	int* zs0 = &grid0[54];
	memcpy(boxdb3, &boxd[3], sizeof boxdb3);
	memcpy(rowdb3, &rowd[3], sizeof rowdb3);
	// now loop over the 24 cells not yet assigned in the band to fill the band use relative cell 
	int ii = -1, free[24];
	uint32_t d;
nextii:
	ii++;
	{
		crcb = tgen_band_cat[ii];//cell_row_col_box one of the 24 cells to fill
		register int fr = cd[crcb[2]] & bd[crcb[3]] & rd[crcb[1]];
		if (!fr)goto back;
		free[ii] = fr;
	}
	goto next_first;

next:// erase previous fill and look for next
	crcb = tgen_band_cat[ii];
	d = zs0[crcb[0]];
	bit = 1 << d;
	rd[crcb[1]] ^= bit; cd[crcb[2]] ^= bit; bd[crcb[3]] ^= bit;
	if (!free[ii])goto back;
	{
	next_first:
		crcb = tgen_band_cat[ii];// be sure to have the good one
		bitscanforward(d, free[ii]);
		bit = 1 << d;
		free[ii] ^= bit;
		zs[crcb[0] + 54] = (char)(d + '1');
		zs0[crcb[0]] = d;
		rd[crcb[1]] ^= bit; cd[crcb[2]] ^= bit; bd[crcb[3]] ^= bit;
		if (ii < 23) goto nextii;
		// this is a valid band, check if canonical 
		int ir = bandminlex.Getmin(&grid0[54], &pband3, 0);
		if (ir < 0) return; //would be bug  never seen
		it16_3 = pband3.i416;	i3t16 = t416_to_n6[it16_3];
		if (op.bx3 < 416 && (op.bx3 != i3t16)) goto next;
		int locdiag = 0;
		//if (nb12 == 469978)locdiag = 2;
		//if (nb12 == 9365+1)locdiag = 2;
			//fulldiag(0);
		if (op.p1) {//Bx c<=a<=b
			if (i3t16 > i1t16) goto next;
			if (sgo.bfx[6]) {
				if(i3t16 != i2t16)goto next;
			}
			if (locdiag) {
				for (int i = 0; i < 27; i++) cout << zs0[i] + 1;
				cout << "new b3 i3t16 " << i3t16 << endl;
			}
			irtw2 = tww1.CheckDiagP1(*this);
			//if (locdiag) cout << "back irtw2 =" << irtw2  << endl;
			if (irtw2 == 1 )goto next;
			if (F3B_See_NED3()&op.test_ned)
				bands3[nband3++].InitBand3(it16_3, &zsol[54], pband3);
			goto next;
		}
		else if (op.p2b) {//Bx a<=c<=b
			if (i3t16 < i1t16) goto next;
			if (i2t16 < i3t16) goto next;
			irtw2 = tww1.CheckDiagP2b(*this);
			if (irtw2 == 1)goto next;

			goto next;
		}
		else {// op.p2a Bx a<=b<=c
			if (i3t16 <  i2t16) goto next;
			irtw2 = tww1.CheckDiagP2a();
			if (irtw2 == 1)goto next;
			if (sgo.bfx[6]) {
				if (i3t16 != i1t16)goto next;
			}
			if(F3B_See_NED1() & op.test_ned)
				bands3[nband3++].InitBand3(it16_3, &zsol[54], pband3);
			goto next;
		}
	}
back:
	if (--ii >= 0) goto next;
	int nn = 0;
	if (nband3) {
		if (op.out_entry) OutEntry();
		else if (m10 == 1)g17b.Start();// call the process for that entry
	}

}


void BandReShape(int* s, int* d, BANDMINLEX::PERM p);
void BandReOrder(int* d);
void GEN_BANDES_12::F3pass2_See() {
	int ir = bandminlex.Getmin(&grid0[54], &pband3, 0);
	if (ir < 0) {//would be bug  did not come in enumeration
		cerr << "gen band 3 invalid return Getmin" << endl;
		return;
	}
	int it16_3 = pband3.i416;
	ib3check = i3t16 = t416_to_n6[it16_3];
	if (op.bx3 < 416)if (op.bx3 != i3t16) return;

	if (i3t16 < i1t16) return;// not canonical

	if (!op.p2b) {// p2a
		if (i3t16 < i2t16)return;// not canonical (must be in this case
		int locdiag = 0;
		if (locdiag) {
			for (int i = 54; i < 81; i++)cout << grid0[i] + 1;
			cout << " band3 to check i3t16=" << i3t16 << endl;
		}
		pcheck3 = pband3;
		memcpy(&gcheck[54], &grid0[54], 27 * sizeof gcheck[0]);
		//if (Band3Check())return;
		if (locdiag) cout << " valid" << endl; 
	}
	else {// p2b exchanging band 2 band 3
		int locdiag = 0;
		if (locdiag) {
			for (int i = 54; i < 81; i++)cout << grid0[i] + 1;
			cout << " band3 to check i3t16=" << i3t16 << endl;
		}
		if (i3t16 > i2t16)return;
		memcpy(&gcheck[27], &grid0[54], 27 * sizeof gcheck[0]);
		pcheck3 = pband2;
		pcheck2 = pband3;
		ib2check = i3t16;
		ib3check = i2t16;
		if (Band2Check())return;// band 3 must be a valid "band2"
		//if (Band3Check())return;// then band 2 a valid band3
	}
	bands3[nband3++].InitBand3(it16_3, &zsol[54], pband3);
}
inline void GEN_BANDES_12::F3B_See_18() {// one NED return 1 if equal not loaded
	pcheck2 = pband3;	pcheck3 = pband2;
	ib1check = i1t16;	ib2check = i3t16;	ib3check = i2t16;
	ibasecheck = it16;
	memcpy(&gcheck[27], &grid0[54], 27 * sizeof gcheck[0]);
	memcpy(&gcheck[54], &grid0[27], 27 * sizeof gcheck[0]);
	if (Band2Check()) return;
	//if (Band3Check()) return;
	bands3[nband3++].InitBand3(it16_3, &zsol[54], pband3);
}

void GEN_BANDES_12::fulldiag(int mode) {
	for (int i = 0; i < 81; i++) cout << grid0[i] + 1;
	cout << " " << i2t16 << " " << i3t16 << "  irtw2=" << irtw2
		<< " in diag  wb2b1.ntab1b2=" << wb2b1.ntab1b2
		<< "  wb1b2.ntab1b2=" << wb1b2.ntab1b2
		<< " slice "<< (nb12 >> 6) << " mode="<<mode << endl;
	if (!mode)return;
	tww1.Dump(" diag status");

	for (int i = 0; i < 6; i++) {
		tww2.PushP1ToMinimal(*this, grid0, i);// perm 
		for (int i = 0; i < 81; i++) cout << tww2.zs0[i] + 1;
		cout << " minimal perm " << i << endl;
		//if (tww.Compare(tww2.zs0) < 0) return;
	}
	if (mode<2)return;
	cout <<"min on diag"<<endl;
	for (int i = 0; i < 6; i++) {
		tww2.PushP1ToMinimal(*this, tww1.zs0, i);// perm 
		for (int i = 0; i < 81; i++) cout << tww2.zs0[i] + 1;
		cout << " minimal perm " << i << endl;
		if (tww.Compare(tww2.zs0) < 0) return;
	}

}
int GEN_BANDES_12::TWW1::CheckDiagP2a() {
	GEN_BANDES_12& o = genb12;
	InitD(o.grid0);
	if (ibx[0] < o.i1t16) return 1;	if (ibx[0] > o.i1t16) return 0;
	if (ibx[1] < o.i2t16) return 1;	if (ibx[1] > o.i2t16) return 0;
	if (ibx[2] < o.i3t16) return 1;	if (ibx[2] > o.i3t16) return 0;
	// same both ways 
	MorphToB1First();
	int ir = Compare(o.grid0);
	if (ir > 0) return 1; // diag smaller
	return 2;
}
int GEN_BANDES_12::TWW1::CheckDiagP2aMore_abc() {
	GEN_BANDES_12& o = genb12;
	// same both ways 3 different Bx
	if (!o.n_auto_b1) return 2;
	for (int i = 0; i < o.n_auto_b1; i++) {
		perm_ret = o.t_auto_b1[i];
		MorphToB1();
		if(Compare(o.grid0)>0)return 0; // diag smaller
	}
	return 2;
}
int GEN_BANDES_12::TWW1::CheckDiagP1More_abc() {
	GEN_BANDES_12& o = genb12;
	// same both ways 3 different Bx
	if (!o.n_auto_b1) return 2;
	for (int i = 0; i < o.n_auto_b1; i++) {
		perm_ret = o.t_auto_b1[i];
		MorphToB1();
		if (Compare(o.grid0) > 0)return 0; // diag smaller
	}
	return 2;
}


int GEN_BANDES_12::TWW::IsBelowP2a(){
	GEN_BANDES_12& o = genb12;
	int* zz = &o.grid0[54], * z0 = &zs0[54];
	//if (o.wb1b2.ntab1b2)	if (o.wb1b2.Check3(zz, z0)) return 1;
	if (o.wb2b1.ntab1b2) { //first  morph b3 to min b2b1
		BANDMINLEX::PERM& p = o.pband2;
		int zt[27];		BandReShape(z0, zt, p);
		if (o.wb2b1.Check3(zz, zt)) return 1;
	}
	if (o.i3t16 != o.i2t16) return 0;
	for (int imorph = 0; imorph < o.n_auto_b1; imorph++) {
		int* z = &zs0[54];
		BANDMINLEX::PERM& p = o.t_auto_b1[imorph];
		SKT_MORPHTOP
			int ir = G17ComparedOrderedBand(&o.grid0[27], band);
		if (ir > 1) continue;		if (ir == 1) return 1;
		// same b2 must check b3
		memcpy(z, &zs0[54], sizeof z);
		{	int* z = &zs0[27];
		SKT_MORPHTOP
			if (G17ComparedOrderedBand(&o.grid0[54], band) == 1) return 1;
		}
	}
	return 0;
}

int GEN_BANDES_12::B1B2::Diag3(int* z) {// z is band3 to check
	for (int imorph = 0; imorph < ntab1b2; imorph++) {
		BANDMINLEX::PERM& p = tab1b2[imorph];
		SKT_MORPHTOP
			BandReOrder(band);
		for (int i = 0; i < 27; i++) cout << band[i] + 1;
		cout << " automorph b1b2 i=" << imorph << endl;
	}
	return 0;
}
int GEN_BANDES_12::fulldiagbug() {
              //214865937768932415935714268341598672576241893892673541;
	char* zt = "214968573368574912975213468531642897742891635896735241";
	for (int i = 27; i < 81; i++) 	if ((zt[i - 27]-'1') != grid0[i]) return 0;

	for (int i = 0; i < 27; i++) cout << grid0[i] + 1; cout << endl;
	for (int i = 27; i < 54; i++) cout << grid0[i] + 1; cout << endl;
	for (int i = 54; i < 81; i++) cout << grid0[i] + 1;
	cout << " " << i2t16 << " " << i3t16 << " in diag  b2b1.n=" << wb2b1.ntab1b2
		<< "  b1b2.n=" << wb1b2.ntab1b2 << " nb12 " << nb12 << " slice " << (nb12 >> 6)
		<< endl << "this is the current status" << endl;;
	n_auto_p1 = tblnauto[it16_2]; //ia 0-415 not index
	t_auto_p1 = &automorphsp[tblnautostart[it16_2]];

/*
		for (int imorph = 0; imorph < n_auto_p1; imorph++) {
			BANDMINLEX::PERM& p = t_auto_p1[imorph];
			int* z = &grid1[27];// morph the band
			SKT_MORPHTOP
				register int ir = G17ComparedOrderedBand(&grid1[27], band);
			if (ir > 1) continue;
			if (ir == 1) {// new mini
				BandReOrder(band);
				memcpy(&zs1[27], band, sizeof band);
				BandReShape(&grid1[54], &zs1[54], p);
				continue;
			}
			int zt[27];// now possible new min band3
			BandReShape(&grid1[54], zt, p);
			ir = BandCompare(zt, &zs1[54]);
			if(ir<1)memcpy(&zs1[54], zt, sizeof zt);
		}
		if (locdiag) GridDump(zs1," P1 min");
*/

	tww2.PushP1ToMinimal(*this, grid0, 3);// perm bca
	tww2.DumpP1("expected min");
	n_auto_p1 = n_auto_b1; //ia 0-415 not index
	t_auto_p1 = t_auto_b1;
	wb1b2.ntab1b2 = wb2b1.ntab1b2 = 0;
	if (n_auto_p1) {// rebuild aoto morphisms
		wb1b2.Init(grid1);
		cout << " go check p1" << endl;
		int* z = &grid1[27];
		for (int imorph = 0; imorph < n_auto_p1; imorph++) {
			BANDMINLEX::PERM& p = t_auto_p1[imorph];
			SKT_MORPHTOP
				int ir = G17ComparedOrderedBand(&grid1[27], band);
			if (!ir ) wb1b2.tab1b2[wb1b2.ntab1b2++] = p;
		}
		cout << " wb1b2.ntab1b2="<< wb1b2.ntab1b2 << endl;
	}
	if (wb1b2.ntab1b2) {
		cout << "try automorph b1b2" << endl;
		wb1b2.Diag3(&tww2.zs0[54]);
	}
	return 1;
}
int GEN_BANDES_12::fulldiagbugP1() {
	char* zt = "274915863365874912918263475546738291732591648891642537";
	//          231764895765298431894531267312945678578613942946872513";
	for (int i = 27; i < 81; i++) 	if ((zt[i - 27] - '1') != grid0[i]) return 0;
	cout << "NED3" << endl;
	for (int i = 0; i < 27; i++) cout << grid0[i] + 1; cout << endl;
	for (int i = 27; i < 54; i++) cout << grid0[i] + 1; cout << endl;
	for (int i = 54; i < 81; i++) cout << grid0[i] + 1;
	cout << " " << i2t16 << " " << i3t16 << " in diag  nb12 " << nb12 << " slice " << (nb12 >> 6) << endl;
	return 1;
}

int GEN_BANDES_12::F3B_See_NED1() {// pass 2a a<=b<=c
	int locdiag = 0;
	if (i3t16 == i2t16 && grid0[27] != 1) return 0;
	if (wb1b2.ntab1b2 && wb1b2.Check3(&grid0[54], &grid0[54])) return 0;
	//if (sgo.bfx[6])fulldiag(1);
	if (i1t16 < i2t16 && i2t16 < i3t16) {
		if (irtw2 != 2) return 2;
		return tww1.CheckDiagP2aMore_abc();
	}
	if (i1t16 == i3t16) {//	three bands = 
		if (sgo.bfx[6] == 2) 	if (!fulldiagbug()) return 0;
		else locdiag = 1;
		for (int i = 0; i <= 5; i++) {
			if (locdiag)cout << "try morep2a i=" << i << endl;
			if (tww2.IsBelowMoreP2a(*this, grid0, i,locdiag)) {
				if (locdiag)cout << "exit ret code="
					<< tww2.IsBelowMoreP2a(*this, grid0, i)	<< endl;
				return 0;
			}
			if (locdiag)cout << "ok morep2a i=" << i << endl;
		}
		if (irtw2 != 2) return 4;
		for (int i = 0; i <= 5; i++)
			if (tww2.IsBelowMoreP2a(*this, tww1.zs0, i))return 0;
		return 4;
	}
	tww.Init(grid0);
	if (tww.IsBelowP2a()) return 0;
	if (irtw2 == 2) {// see diagonal
		if (tww2.IsBelowMoreP2a(*this, tww1.zs0, 5))return 0;
		if (i1t16 == i2t16 &&
			tww2.IsBelowMoreP2a(*this, tww1.zs0, 4))return 0;
		if (i3t16 == i2t16 &&
			tww2.IsBelowMoreP2a(*this, tww1.zs0, 3))return 0;
	}
	if (i1t16 == i2t16) return 8;
	return 1;
}

int GEN_BANDES_12::F3B_See_NED2() {// // pass 2b a<=c<=b 
	if (sgo.bfx[6])fulldiag(2);
	if (i2t16 == i3t16) {// same 3 bands, more perms to see 
		// first get minimal can be abc=> cab or cba  or acb or bca
		int grid1[81], mint = 0; grid1[0] = 10;
		for (int i = 0; i < 2; i++) {
			tww2.PushP1ToMinimal(*this, grid0, i);// perm 
			int ir = tww2.Compare0(grid1);
			if (ir > 0) {
				mint = 1 << i;
				memcpy(grid1, tww2.zs0, sizeof grid1);
			}
			else if (ir == 0)mint |= 1 << i;
		}
		for (int i = 0; i < 81; i++) {
			if (grid0[i] > grid1[i]) break;
			if (grid0[i] < grid1[i]) return 0;
		}
		for (int i = 2; i <= 4; i++) {
			tww2.PushP1ToMinimal(*this, grid0, i);// perm 
			if (tww2.Compare0(grid1) > 0)  return 0;
		}
		if (irtw2)for (int i = 0; i < 6; i++) {
			tww2.PushP1ToMinimal(*this, tww1.zs0, i);// perm 
			if (tww2.Compare0(grid1) > 0)  return 0;
		}
		if (mint == 2 && wb2b1.ntab1b2) return 0;
		return 4;
	}
	if (wb1b2.ntab1b2 && wb1b2.Check3(&grid0[54], &grid0[54])) return 0;
	if (i3t16 < i1t16 && i1t16 < i2t16) {// bands 1+2 checked minimal
		if (!irtw2) return 2;
		tww.PushP1ToMinimal(*this, grid0);// expected min
		tww2.PushP1ToMinimal(*this, tww1.zs0, 0);//no  perm 
		if (tww.Compare(tww2.zs0) < 0) return 0;
		return 2;
	}
	if (i1t16 == i3t16) {// i3t16;i1t16 must be minimal
		// if auto morphs band1 band2 keep only lower band 3
		// i2t16 is highest min  can be 0/cab or 2/acb 
		tww.PushP1ToMinimal(*this, grid0);// expected min
		//fulldiag(1);
		tww2.PushP1ToMinimal(*this, grid0, 2);// perm  acb
		if (tww.Compare(tww2.zs0) < 0) {
			if (!wb2b1.ntab1b2)	return 0;
			memcpy(tww.zs0, tww2.zs0, sizeof tww.zs0);
		}
		//cout << " ok to go" << endl;
		if (irtw2 == 2) {	// tww1 is in pass 1 mode
			//cout << "i1t16 == i3t16 and irtw2 == 2" << endl;
			tww2.PushP1ToMinimal(*this, tww1.zs0, 2);// perm  acb
			if (tww.Compare(tww2.zs0) < 0) return 0;
			tww2.PushP1ToMinimal(*this, tww1.zs0);// perm  cab
			if (tww.Compare(tww2.zs0) < 0) return 0;
		}
		return 8;
	}
	tww.Init(grid0); // check auto morphs direct 
	if (tww.IsBelowP1(*this)) return 0;


	tww.InitPermb1b2(grid0); //now i2t16 = i1t16
	if (tww.IsBelowP1(*this)) return 0;
	if (irtw2 == 2) {	// tww1 is in pass 1 mode 
		if (tww1.IsBelowP1(*this)) return 0;
		tww1.PermB(0, 1);		tww1.MorphToB1First();
		int ir = tww1.Compare(grid0);
		if (ir > 0) return 0; // diag smaller
	}
	return 1;
}

int GEN_BANDES_12::NED3_Build_Check_Morphs(int locdiag) {
	// clear band3 redundancy
	if (wb1b2.ntab1b2 && wb1b2.Check3(&grid0[54], &grid0[54])) return 0;	// rebuild nauto b1b2 b2b1 on band1 auto morphs 
	tww.MorphToP1(*this);
	n_auto_p1 = tblnauto[it16_3]; //ia 0-415 not index
	t_auto_p1 = &automorphsp[tblnautostart[it16_3]];
	if (n_auto_p1) {// rebuild aoto morphisms
		wb1b2.Init(grid1);
		int zs1[81]; memcpy(zs1, grid1, sizeof zs1);
		// move to mini
		for (int imorph = 0; imorph < n_auto_p1; imorph++) {
			BANDMINLEX::PERM& p = t_auto_p1[imorph];
			int* z = &grid1[27];// morph the band
			SKT_MORPHTOP
				BandReOrder(band);
			if (locdiag) BandDump(band, " morph b2 of min");
			register int ir = G17ComparedOrderedBand(&zs1[27], band);
			if (ir > 1) continue;
			if (ir == 1) {// new mini 
				memcpy(&zs1[27], band, sizeof band);
				BandReShape(&grid1[54], &zs1[54], p);
				continue;
			}
			int zt[27];// now possible new min band3
			BandReShape(&grid1[54], zt, p);
			ir = BandCompare(zt, &zs1[54]);
			if (ir < 1) 	memcpy(&zs1[54], zt, sizeof zt); 
			
		}
		if (locdiag) GridDump(zs1," P1 min");
		memcpy(grid1, zs1, sizeof grid1);
	}
	if (i2t16 == i1t16) {// grid1 Bx2=Bx3
		if (grid1[27] != 1) return 0;
		if (tww.IsBelowGrid1(*this, grid1, 0))return 0; // acb
	}
	if (i3t16 == i1t16) {	// grid1 Bx1=Bx2
		if (tww.IsBelowGrid1(*this, grid1, 1))return 0;   // bac
	}
	return 1;
}


int GEN_BANDES_12::F3B_See_NED3() {//  pass 1 c<=a<=b 
	int locdiag = 0;
	//if (wb2b1.ntab1b2) locdiag = 2;
	//if ((nb12 >> 6) == 35) locdiag = 2;
	if (sgo.bfx[6]) {
		if (fulldiagbugP1())locdiag = 2;
		//else return 0;
		if (locdiag) fulldiag(0);
	}
	if (!NED3_Build_Check_Morphs(locdiag)) return 0;

	if (i3t16 < i1t16 && i1t16 < i2t16) {
		//if (!irtw2) return 2;
		if (irtw2 != 2) return 2;
		return tww1.CheckDiagP1More_abc();
	}

	if (i2t16 == i3t16) {// same 3 bands, more perms to see 
		if (locdiag) cout << "nhr" << endl;
		for (int i = 0; i <= 5; i++)
			if (tww.IsBelowGrid1(*this, grid1, i))return 0;
		if (irtw2) cout << "diag to check" << endl;
		if (irtw2)for (int i = 0; i < 6; i++) {
			if (tww.IsBelowGrid1(*this, grid1, i))return 0;
		}

		for (int i = 27; i < 81; i++) fout1 << grid1[i] + 1;
		fout1 << ";";
		for (int i = 27; i < 81; i++) fout1 << grid0[i] + 1;
		fout1 << " nb12="<< nb12  << endl;
		return 0;
		return 4;
	}
	if (i1t16 == i3t16) {// i3t16;i1t16 must be minimal
		// i2t16 is highest min  can be 0/cab or 2/acb 
		//tww.PushP1ToMinimal(*this, grid0);// expected min
		//fulldiag(1);
		//tww2.PushP1ToMinimal(*this, grid0, 2);// perm  acb
		//if (tww.Compare(tww2.zs0) < 0) {
			//if(!wb2b1.ntab1b2)	return 0;
			//memcpy(tww.zs0, tww2.zs0, sizeof tww.zs0);
		//}
		//cout << " ok to go" << endl;
		if (irtw2 == 2) {	// tww1 is in pass 1 mode
			if (sgo.bfx[6])fulldiag(5);
			//cout << "i1t16 == i3t16 and irtw2 == 2" << endl;
			tww2.PushP1ToMinimal(*this, tww1.zs0, 2);// perm  acb
			if (tww.Compare(tww2.zs0) < 0) return 0;
			tww2.PushP1ToMinimal(*this, tww1.zs0);// perm  cab
			if (tww.Compare(tww2.zs0) < 0) return 0;
		}
		return 8;
	}
	tww.Init(grid0); // check auto morphs direct 
	if (tww.IsBelowP1(*this)) return 0;


	tww.InitPermb1b2(grid0); //now i2t16 = i1t16
	if (tww.IsBelowP1(*this)) return 0;
	if (irtw2 == 2) {	// tww1 is in pass 1 mode 
		if (tww1.IsBelowP1(*this)) return 0;
		tww1.PermB(0, 1);		tww1.MorphToB1First();
		int ir = tww1.Compare(grid0);
		if (ir > 0) return 0; // diag smaller
	}
	return 1;
}


 