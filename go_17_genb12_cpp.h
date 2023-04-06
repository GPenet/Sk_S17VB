

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
	myband1.InitG12(it16);

	memcpy(grid0, myband1.band0, sizeof myband1.band0);
	memcpy(gcheck, myband1.band0, sizeof myband1.band0);

	strcpy(zsol, myband1.band);
	n_auto_b1 = bandminlex.GetAutoMorphs(it16, t_auto_b1);
	for (int i = 0; i < 9; i++) // init columns status
		cold[i] = 0x1ff ^ myband1.gangster[i];
	memcpy(coldf, cold, sizeof coldf);
	zsol[27] = 0;
	if(op.ton)
	cout << myband1.band << "i1t16=" << i1t16 << " it16=" << it16
		<< " n auto morphs=" << n_auto_b1 << endl;
	if (n_auto_b1) {
		int * zs0 = grid0;
		for (int imorph = 0; imorph < n_auto_b1; imorph++) {
			BANDMINLEX::PERM &p = t_auto_b1[imorph];
			int band[27];// morph the band
			for (int i = 0; i < 9; i++) {
				band[i] = p.map[zs0[p.cols[i]]];
				band[i + 9] = p.map[zs0[p.cols[i] + 9]];
				band[i + 18] = p.map[zs0[p.cols[i] + 18]];
			}
		}
	}
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
int GEN_BANDES_12::Band3Check() {
	if (i1t16 == ib3check && ib3check == ib2check) {// 3 bands equal use diagonal test 
		BANDMINLEX::PERM* p = minlexusingbands.pout;
		p[0].InitBase(i1t16);
		p[1] = pcheck2;
		p[2] = pcheck3;
		if (minlexusingbands.IsLexMinDirect(gcheck, i1t16, t_auto_b1, n_auto_b1))
			return 1;
		return 0;
	}
	{
		//========================== morphs on b1b2 base test
		if (n_auto_b1b2) {// still direct automorphism b1b2
			for (int imorph = 0; imorph < n_auto_b1b2; imorph++) {
				BANDMINLEX::PERM &p = t_auto_b1b2[imorph];
				int b23[3][9];
				// direct
				for (int i = 0; i < 3; i++) {// band 3 only
					register int * rrd = b23[i], *rro = &gcheck[54 + 9 * i];
					for (int j = 0; j < 9; j++)		rrd[j] = p.map[rro[p.cols[j]]];
				}
				if (G17ComparedOrderedBand(&gcheck[54], b23[0]) == 1)
					return 1;
			}
		}
		//=========================== perm b1b2 and base test (b1=b2)
		if (n_auto_b2b1) {// possible lower band3 with a perm band1 band2
			int b23[3][9];//first morph to band 2 min lexical
			for (int i = 0; i < 3; i++) {// rows 4 to 9 as of band 2 perm
				register int * rrd = b23[i], *rro = &gcheck[54 + 9 * i];
				for (int j = 0; j < 9; j++)
					rrd[j] = pcheck2.map[rro[pcheck2.cols[j]]];
			}
			for (int imorph = 0; imorph < n_auto_b2b1; imorph++) {// then apply auto morphs 
				BANDMINLEX::PERM &pp = t_auto_b2b1[imorph];
				int b23_a[3][9];
				for (int i = 0; i < 3; i++) {
					register int * rrd = b23_a[i], *rro = b23[i];
					for (int j = 0; j < 9; j++)		rrd[j] = pp.map[rro[pp.cols[j]]];
				}
				if (G17ComparedOrderedBand(&gcheck[54], b23_a[0]) == 1)
					return 1;
			}
		}
		//========================= (b2=b3)#b1  perm b2b3 to consider (direct done)
		if (ib3check == ib2check) {// check b3b2 on  auto morphs b1
			if (gcheck[27] - 1) return 1; // must be '2' in r4c1
			for (int imorph = 0; imorph < n_auto_b1; imorph++) {
				BANDMINLEX::PERM &p = t_auto_b1[imorph];
				int b23[6][9];
				for (int i = 0; i < 6; i++) {// rows 4 to 9 from band 2
					register int * rrd = b23[i], *rro = &gcheck[27 + 9 * i];
					for (int j = 0; j < 9; j++)		rrd[j] = p.map[rro[p.cols[j]]];
				}
				int ir = G17ComparedOrderedBand(&gcheck[27], b23[3]);
				if (ir == 1)return 1;
				if (ir < 1 && G17ComparedOrderedBand(&gcheck[54], b23[0]) == 1)return 1;
			}
		}
		//============================ 
		if (minlexusingbands.IsLexMinDiagB(gcheck, i1t16, ib2check, ib3check, t_auto_b1, n_auto_b1))
		   return 1;
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
		pcheck2 = pband2;
		it16_2 = pband2.i416;
		ib2check = i2t16 = t416_to_n6[it16_2];
		if (i2t16 < i1t16)goto next;// not canonical
		if (op.p2b||op.p1)memcpy(&gcheck[54], zs0, 27 * sizeof gcheck[0]);
		else {
			memcpy(&gcheck[27], zs0, 27 * sizeof gcheck[0]);
			if (Band2Check())goto next;// do nothing if p2b or p1
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
int GEN_BANDES_12::ValidBand2() {
	if (g17b.aigstop)return 1;
	myband2.InitBand2_3(it16_2, &zsol[27], pband2);
	//_______________________ std process
	if (modeb12 < 10) {
		nband3 = 0;
		if((nb12 >> 6) < op.first) return 0;// here restart value, kept untouched if no band 3 found
		{// print a restart point every 64 bands 1+2 seen
			uint64_t w = genb12.nb12, w1 = w >> 6;
			w &= 63;
			if (w == 0) {
				long tfin = GetTimeMillis();
				cout << "next slice to use=\t" << w1 << "\tmil=" << (tfin - sgo.tdeb) / 1000 << "\tnb2=" << p_cpt2g[0] << endl;
			}
		}
		if (((nb12 >> 6) > op.last)) return 1;
		ValidInitGang();// also called from existing 17 in test
		if(F17Novalid1_2())return ((nb12 >> 6) > op.last);
		if (op.p1)Find_band3B_pass1B();
		else Find_band3B();
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
	
		//if (op.p1)Find_band3B_pass1(0);
		if (op.p1)Find_band3B_pass1B(0);
		else Find_band3B(0);
		if (nband3) {	p_cpt[0]++;	p_cpt[1] += nband3;	}
		if (((nb12 >> 6) > op.last)) return 1;
	}
	return 0;
}
void GEN_BANDES_12::OutEntry() {
	p_cpt2g[1] += nband3;   // update bands3 count
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
	//BANDMINLEX::PERM pout;
	register int  *crcb, bit;
	nband3 =0;
	int *rd = rowdb3, *cd = cold, *bd = boxdb3; // to updates rows cols boxes
	char * zs = zsol;
	int * zs0 = &grid0[54];
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
		int ir = bandminlex.Getmin(zs0, &pband3, 0);
		if (ir < 0) {//would be bug  did not come in enumeration
			cerr << "gen band 3 invalid return Getmin" << endl;
			return;
		}	
		int it16_3 = pband3.i416;
		ib3check=i3t16 = t416_to_n6[it16_3];
		if (op.bx3 < 416)if (op.bx3 != i3t16) goto next;

		if (i3t16 < i1t16)goto next;// not canonical
		//if (op.b2 && op.ton) {
			//for (int i = 0; i < 27; i++)cout << zs0[i] + 1;
			//cout <<"seen" <<i2t16<<" "<<i3t16<<endl;
		//}
		if (!op.p2b) {// p2a
			if (i3t16 < i2t16)goto next;// not canonical (must be in this case
			pcheck3 = pband3;
			memcpy(&gcheck[54], zs0, 27 * sizeof gcheck[0]);
			if (Band3Check())goto next;
		}
		else {// p2b exchanging band 2 band 3
			if (i3t16 > i2t16)goto next;
			memcpy(&gcheck[27], zs0, 27 * sizeof gcheck[0]);
			pcheck3 = pband2;
			pcheck2 = pband3;
			ib2check = i3t16;
			ib3check = i2t16;
			if (Band2Check())goto next;// band 3 must be a valid "band2"
			if (Band3Check())goto next;// then band 2 a valid band3
		}
		bands3[nband3++].InitBand3(it16_3, &zs[54], pband3);
		goto next;
	}
back:
	if (--ii >= 0) goto next;
	if (m10 != 1)return;
	if (nband3) {
		if (op.out_entry)OutEntry();
		else g17b.Start();// call the process for that entry
	}
}
void GEN_BANDES_12::Find_band3B_pass1B(int m10) {
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
		F3B_See();
		goto next;
	}
back:
	if (--ii >= 0) goto next;
	//if (m10 != 1)return;
	if (nband3) {
		if (op.out_entry) OutEntry();
		else if (m10 == 1)g17b.Start();// call the process for that entry
	}

}

void BandReShape(int* s, int* d, BANDMINLEX::PERM p) {
	int* pc = p.cols;
	//uint8_t* pc = p.cols;
	for (int irow = 0; irow < 3; irow++) {
		int drow = 9 * irow;
		for (int icol = 0; icol < 9; icol++)
			d[drow + icol] = p.map[s[drow + pc[icol]]];
	}
	int temp[9];// now reorder 
	if (d[0] > d[9]) {
		memcpy(temp, &d[0], sizeof temp);
		memcpy(&d[0], &d[9], sizeof temp);
		memcpy(&d[9], temp, sizeof temp);
	}
	if (d[0] > d[18]) {
		memcpy(temp, &d[0], sizeof temp);
		memcpy(&d[0], &d[18], sizeof temp);
		memcpy(&d[18], temp, sizeof temp);
	}
	if (d[9] > d[18]) {
		memcpy(temp, &d[9], sizeof temp);
		memcpy(&d[9], &d[18], sizeof temp);
		memcpy(&d[18], temp, sizeof temp);
	}
}
void BandReOrder(int* d) {
	int temp[9];// now reorder 
	if (d[0] > d[9]) {
		memcpy(temp, &d[0], sizeof temp);
		memcpy(&d[0], &d[9], sizeof temp);
		memcpy(&d[9], temp, sizeof temp);
	}
	if (d[0] > d[18]) {
		memcpy(temp, &d[0], sizeof temp);
		memcpy(&d[0], &d[18], sizeof temp);
		memcpy(&d[18], temp, sizeof temp);
	}
	if (d[9] > d[18]) {
		memcpy(temp, &d[9], sizeof temp);
		memcpy(&d[9], &d[18], sizeof temp);
		memcpy(&d[18], temp, sizeof temp);
	}
}


//============== pass1 filters
void GEN_BANDES_12::F3B_See() {
	int ir = bandminlex.Getmin(&grid0[54], &pband3, 0);
	if (ir < 0)return; // bug  never seen
	it16_3 = pband3.i416;
	i3t16 = t416_to_n6[it16_3];
	if (op.bx3 < 416)if (op.bx3 != i3t16) return;
	if (i3t16 > i2t16) return;// direct not a pass 1 
	// reverse case in p2b one NED in even p1 mode
	if (op.t18) if (i3t16 >= i1t16)F3B_See_18();
	// always b3 <= b1 one NED
	if (i3t16 <= i1t16)  F3B_See_Com();
}

inline void GEN_BANDES_12::F3B_See_18() {// one NED return 1 if equal not loaded
	pcheck2 = pband3;	pcheck3 = pband2;
	ib1check = i1t16;	ib2check = i3t16;	ib3check = i2t16;
	ibasecheck = it16;
	memcpy(&gcheck[27], &grid0[54], 27 * sizeof gcheck[0]);
	memcpy(&gcheck[54], &grid0[27], 27 * sizeof gcheck[0]);
	if (Band2Check()) return;
	if (Band3Check()) return;
	bands3[nband3++].InitBand3(it16_3, &zsol[54], pband3);
}
void GEN_BANDES_12::F3B_See_Com() {// one NED return 1 if equal not loaded
	if (n_auto_b1) { 	//redundancy if not b1b2 minimal
		int* z = &grid0[27];
		for (int imorph = 0; imorph < n_auto_b1; imorph++) {
			BANDMINLEX::PERM& p = t_auto_b1[imorph]; 	SKT_MORPHTOP
				int ir = G17ComparedOrderedBand(z, band);
			if (ir == 1)	return;
		}
	}
	char wb1[28];// band in char mode
	int wb1i[27]; // band in 0-8 integer mode
	strcpy(wb1, "12345678945");
	strncpy(&wb1[11], t416[it16_3], 16);
	wb1[27] = 0;
	for (int i = 0; i < 27; i++) wb1i[i] = wb1[i] - '1';
	memcpy(gw, wb1i, sizeof wb1i);
	BandReShape(grid0, &gw[27], pband3);
	BandReOrder(&gw[27]);
	BandReShape(&grid0[27], &gw[54], pband3);
	BandReOrder(&gw[54]);
	// push it to minimal b1b2
	ib1check = i3t16;	ib2check = i1t16;	ib3check = i2t16;
	ibasecheck = it16_3;
	F3B_See_Com_GetCFX();

}

int GEN_BANDES_12::F3B_See_Com_FilterDiag(){
	//must do here a stack control
	int  zs0d[81];
	for (int i = 0; i < 81; i++) {
		zs0d[i] = gw[C_transpose_d[i]];
	}
	BANDMINLEX::PERM perm_ret;
	bandminlex.Getmin(zs0d, &perm_ret);
	int ib1d = perm_ret.i416, ib1ad = t416_to_n6[ib1d];
	if (ib1ad < i3t16) return 0;
	bandminlex.Getmin(&zs0d[27], &perm_ret);
	int ib2d = perm_ret.i416, ib2ad = t416_to_n6[ib2d];
	if (ib2ad < i3t16) return 0;
	bandminlex.Getmin(&zs0d[54], &perm_ret);
	int ib3d = perm_ret.i416, ib3ad = t416_to_n6[ib3d];
	if (ib3ad < i3t16) return 0;
	if (ib1ad == i3t16) {
		if (ib2ad < i1t16 || ib3ad < i1t16) return 0;
	}
	else if (ib2ad == i3t16) {
		if (ib1ad < i1t16 || ib3ad < i1t16) return 0;
	}
	else if (ib3ad == i3t16) {
		if (ib1ad < i1t16 || ib2ad < i1t16) return 0;
	}
	return ib1ad;
}

int GEN_BANDES_12::Band2_3CheckNoauto() {
	BANDMINLEX::PERM* t_autom = automorphsp;
	 if (ib1check == ib2check) {// must test the perm
		bandminlex.Getmin(&gw[27], &pcheck2, 0);// re do p2check  
		// morph b1 to b2 see if lower
		int ir,ir2;
		BANDMINLEX::PERM p = pcheck2; 
		{
			int* z = gw;// morph the band
			SKT_MORPHTOP
			ir = G17ComparedOrderedBand(&gw[27], band);
		}
		if (ir == 1) return 0;// don't do the perm if still lower
		if(!ir){
			int* z = &gw[54]; SKT_MORPHTOP
			ir2 = G17ComparedOrderedBand(&gw[54], band);
			if( ir2==1) return 0;
		}
	}
	else if (ib2check == ib3check) {// check perm b2 b3
		if (gw[27] - 1) return 0; //here  must be '2' in r4c1
	}
	if (ib1check == ib3check) {
		 if (gw[27] - 1) return 0; //here  must be '2' in r4c1
		 int ir = G17ComparedOrderedBand(&gw[27], &grid0[27]);
		 if (ir == 1) return 0;		
		 if (!F3B_See_Com_FilterDiag()) return 0;;
		 // try also band3 as first perm 3 1 2 to see if lower
		 bandminlex.Getmin(&gw[54], &pcheck3, 0);// redo perm
		 {
			 int* z = gw;
			 BANDMINLEX::PERM& p = pcheck3; SKT_MORPHTOP
				 int ir = G17ComparedOrderedBand(&gw[27], band);
			 if (ir == 1)  return 0;			 
		 }
		 // try also band3 as first perm 3 2 1 to see if lower
		 {
			 int* z = & gw[27];
			 BANDMINLEX::PERM& p = pcheck3; SKT_MORPHTOP
				 int ir = G17ComparedOrderedBand(&gw[27], band);
			 if (ir == 1) return 0;			
		 }
	 }
	int ir = minlexusingbands.IsLexMinDiagB(gw, ib1check, ib2check, ib3check, t_autom, 0);
	if (ir > 1) return 0;// never seen	
	if (ir)		return 0;
	return 1;// good to process
}

void  GEN_BANDES_12::F3B_See_Com_GetCFX() {
	int na = tblnauto[ibasecheck]; //ia 0-415 not index
	if (!na) { 
		if (Band2_3CheckNoauto()) {// good if no auto morph
			if (op.ton == 3) {
				for (int i = 0; i < 81; i++)fout1 << gw[i] + 1;
				fout1 << ";" << i1t16 << ";" << i2t16 << ";"
					<< i3t16 << " p_cpt2g[91] " << p_cpt2g[91] << endl;
			}
			else 	bands3[nband3++].InitBand3(it16_3, &zsol[54], pband3);
		}
		return;
	}
	int  band2min[27],band3min[27];
	memcpy(band2min, &gw[27], sizeof band2min);
	memcpy(band3min, &gw[54], sizeof band3min);
	BANDMINLEX::PERM* t_autom = &automorphsp[tblnautostart[it16_3]];
	int tmini[108], nmini = 1; tmini[0] = -1;
	for (int imorph = 0; imorph < na; imorph++) {
		int* z = &gw[27];// morph the band
		BANDMINLEX::PERM p = t_autom[imorph]; SKT_MORPHTOP
			int ir = G17ComparedOrderedBand(band2min, band);
		if (ir > 1) continue;
		if (!ir) { tmini[nmini++] = imorph;		continue; }
		// now a lower 
		nmini = 0;
		tmini[nmini++] = imorph;
		BandReOrder(band);
		memcpy(band2min, band, sizeof band2min);
	}
	if (ib2check == ib3check) {// skip if a lower can be seen using band3
		for (int imorph = 0; imorph < na; imorph++) {
			int* z = &gw[54];// morph the band
			BANDMINLEX::PERM p = t_autom[imorph]; SKT_MORPHTOP
				int ir = G17ComparedOrderedBand(band2min, band);
			if (ir == 1) return;			
			if (!ir) {
				int* z = &gw[27];// check perm
				SKT_MORPHTOP
				int ir = G17ComparedOrderedBand(band2min, band);
				if (ir == 1)return;				
			}
		}
	}
	// here, can not change first band
	
	bandminlex.Getmin(band2min, &pcheck2, 0);// re do p2check  
	n_auto_b2b1 = 0;
	if (ib1check == ib2check) { // if lower using b2 as base, not a CFX
		int* z2 = band2min, zbw[27];
		{
			int* z = gw;
			BANDMINLEX::PERM& p = pcheck2; SKT_MORPHTOP
				int ir = G17ComparedOrderedBand(z2, band);
			if (ir == 1) return;			
			if (!ir)t_auto_b2b1[n_auto_b2b1++].InitBase(ib2check);
			memcpy(zbw, band, sizeof zbw);
		}
		for (int imorph = 0; imorph < na; imorph++) {
			int* z = zbw;
			BANDMINLEX::PERM& p = t_autom[imorph];		SKT_MORPHTOP
				int ir = G17ComparedOrderedBand(z2, band);
			if (ir == 1) return;			
			else if (!ir) {// auto morph b1 b2 store it for later
				t_auto_b2b1[n_auto_b2b1++] = p;
			}
		}
	}
	if (tmini[0] >= 0) {// morph b3 to imorph
		BANDMINLEX::PERM& p = t_autom[tmini[0]];
		BandReShape(&gw[54], band3min, p);
		BandReOrder(band3min);
	}
	if (nmini > 1) {// use smallest b3
		int band3minw[27];
		for (int i = 1; i < nmini; i++) {// morph to imorph
			BANDMINLEX::PERM& pw = t_autom[tmini[i]];
			BandReShape(&gw[54], band3minw, pw);
			BandReOrder(band3minw);
			for (int j = 0; j < 27; j++) {
				int ir = band3minw[j] - band3min[j];
				if (ir > 0) break;
				if (ir < 0) {
					memcpy(band3min, band3minw, sizeof band3min);
					break;
				}
			}
		}
	}
	if (ib1check == ib3check) {// if same bx three bands 
		int ir = G17ComparedOrderedBand(band2min, &grid0[27]);
		if (ir == 1) return;
	}
	memcpy(&gw[27], band2min, sizeof band2min);
	memcpy(&gw[54], band3min, sizeof band3min);

	{	//=========================== perm b1b2 and base test (b1=b2)
		if (n_auto_b2b1) {// possible lower band3 with a perm band1 band2
			int zbw[27];
			{  //first morph to band 2 min lexical
				int* z = band3min;
				BANDMINLEX::PERM& p = pcheck2; SKT_MORPHTOP
					memcpy(zbw, band, sizeof zbw);
			}
			for (int imorph = 0; imorph < n_auto_b2b1; imorph++) {// then apply auto morphs
				int* z = zbw;
				BANDMINLEX::PERM& p = t_auto_b2b1[imorph]; SKT_MORPHTOP
					if (G17ComparedOrderedBand(band3min, band) == 1)	return;
			}
		}
	}

	if (!F3B_See_Com_FilterDiag()) return;
	int ir = minlexusingbands.IsLexMinDiagB(gw, ib1check, ib2check, ib3check, t_autom, na);
	if (ir > 1) return;		
	if (ir) return;	
	if (tblnauto[ibasecheck]) {// check if redundant
		for (int ich = 0; ich < nsgchecked; ich++) {
			int* old = sgchecked[ich], aig = 1;
			for (int j = 0; j < 81; j++) {
				if (old[j] != gw[j]) { aig = 0; break; }
			}
			if (aig)	return; //seen  redundant
		}
		{
			if (nsgchecked > p_cpt2g[20]) 	p_cpt2g[20] = nsgchecked;
			int* d = sgchecked[nsgchecked++];
			memcpy(d, gw, sizeof gw);
		}
	}	  
	if (op.ton ==3) {
		for (int i = 0; i < 81; i++)fout1 << gw[i] + 1;
		fout1 << ";" << i1t16 << ";" << i2t16 << ";"
			<< i3t16<<" p_cpt2g[91] "<< p_cpt2g[91] << endl;
	}
	else 	bands3[nband3++].InitBand3(it16_3, &zsol[54], pband3);	
	return;
}
