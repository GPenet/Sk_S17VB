
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
		<< it16		<< " n auto morphs=" << n_auto_b1 
		<<" mode b12="<<modeb12 << endl;

	ntc = 0;
	BitsInTable32(tc, ntc, cold[0]);// first col 6 digits in table
	for (ip20 = 0; ip20 < 20; ip20++) {//0;ip<20 setup initial values for rows columns
		GetStartB2(ip20);
		Find_band2B();
		if (go_back)return;
	}
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
		it16_2 = pband2.i416;
		i2t16 = t416_to_n6[it16_2];
		if (i2t16 < i1t16)goto next;// not canonical
		wb1b2.ntab1b2 = wb2b1.ntab1b2 = 0;
		if (n_auto_b1) {
			wb1b2.Init(grid0);
			if (wb1b2.Check(*this) == 1)goto next;
			if (op.p2 && (!op.p2b )) {// compute an apply auto morphs 
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
				if(op.opcode!=79)
					cout << "new slice =\t" << w1 << "\tmil=" << (tfin - sgo.tdeb) / 1000
					<< "\tnb2=" << p_cpt2g[0] << "\tnb3=" << p_cpt2g[1] << endl;
			}
		}
		if ((nb12 >> 6) > op.last) return 1;
		if (op.b2 < 416 && (t416_to_n6[it16_2] != op.b2)) return 0;
		ValidInitGang();// also called from existing 17 in test
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
		if (op.p1) {//Bx c<=a<=b
			if (i3t16 > i1t16) goto next;
			irtw2 = tww1.CheckDiagP1(*this);
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
			if (F3B_See_NED2() & op.test_ned)
				bands3[nband3++].InitBand3(it16_3, &zsol[54], pband3);

			goto next;
		}
		else {// op.p2a Bx a<=b<=c
			if (i3t16 <  i2t16) goto next;
			irtw2 = tww1.CheckDiagP2a();
			if (irtw2 == 1)goto next;
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

int GEN_BANDES_12::TWW1::CheckDiagP2a() {
	GEN_BANDES_12& o = genb12;
	InitD(o.grid0);
	if (ibx[0] < o.i1t16) return 1;	if (ibx[0] > o.i1t16) return 0;
	if (ibx[1] < o.i2t16) return 1;	if (ibx[1] > o.i2t16) return 0;
	if (ibx[2] < o.i3t16) return 1;	if (ibx[2] > o.i3t16) return 0;
	// same both ways 
	MorphToB1First();
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
int GEN_BANDES_12::TWW::IsBelowP2a(){
	GEN_BANDES_12& o = genb12;
	int* zz = &o.grid0[54], * z0 = &zs0[54];
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
int GEN_BANDES_12::F3B_See_NED1() {// pass 2a a<=b<=c
	if (i3t16 == i2t16 && grid0[27] != 1) return 0;
	if (wb1b2.ntab1b2 && wb1b2.Check3(&grid0[54], &grid0[54])) return 0;

	if (i1t16 < i2t16 && i2t16 < i3t16) {
		if (irtw2 != 2) return 2;
		return tww1.CheckDiagP2aMore_abc();
	}
	if (i1t16 == i3t16) {//	three bands = 
		for (int i = 0; i <= 5; i++) 
			if (tww2.IsBelowMoreP2a(*this, grid0, i)) return 0;
		if (irtw2 != 2) return 4;
		for (int i = 0; i <= 5; i++)
			if (tww2.IsBelowMoreP2a(*this, tww1.zs0, i))return 0;
		return 4;
	}
	tww.Init(grid0);	if (tww.IsBelowP2a()) return 0;
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

int GEN_BANDES_12::NED2_Build_Check_Morphs() {
	if (wb1b2.ntab1b2 && wb1b2.Check3(&grid0[54], &grid0[54])) return 0;	
	tww.MorphToP2b(*this);
	n_auto_p1 = n_auto_b1; 
	t_auto_p1 = t_auto_b1;
	if (n_auto_p1) {// build min canonical Bx
		wb1b2.Init(grid1);
		int zs1[81]; memcpy(zs1, grid1, sizeof zs1);
		// move to mini
		for (int imorph = 0; imorph < n_auto_p1; imorph++) {
			BANDMINLEX::PERM& p = t_auto_p1[imorph];
			int* z = &grid1[27];// morph the band
			SKT_MORPHTOP
				BandReOrder(band);
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
		memcpy(grid1, zs1, sizeof grid1);
	}
	if (i2t16 == i3t16) {// grid1 Bx2=Bx3
		if (grid1[27] != 1) return 0;
		if (tww.IsBelowGrid1(*this, grid1, 0))return 0; // acb
	}
	if (i3t16 == i1t16 	// grid1 Bx1=Bx2
		&& tww.IsBelowGrid1(*this, grid1, 1))return 0;   // bac
	if (irtw2	// check diag direct against min
		&& tww2.IsBelowGrid1(*this, tww1.zs0, 0))return 0;// acb on tww1.zs0
	return 1;
}
int GEN_BANDES_12::F3B_See_NED2() {// // pass 2b a<=c<=b 
	if (!NED2_Build_Check_Morphs()) return 0;
	if (i1t16 < i3t16 && i3t16 < i2t16)
		if (irtw2)  return 2;
	if (i1t16 == i2t16) {// same 3 bands  
		for (int i = 0; i <= 5; i++)
			if (tww.IsBelowGrid1(*this, grid1, i))return 0;
		if (irtw2)for (int i = 0; i < 6; i++) {
			if (tww.IsBelowGrid1(*this, tww1.zs0, i))return 0;
		}
		return 4;
	}
	if (irtw2 == 2) {// see diagonal
		if (i1t16 == i3t16 &&  // Bx1=Bx2 
			tww2.IsBelowGrid1(*this, tww1.zs0, 6))return 0;
		if (i2t16 == i3t16 && // Bx2=Bx3 
			tww2.IsBelowGrid1(*this, tww1.zs0, 5))return 0;
	}
	if (i1t16 == i3t16)		return 8;
	return 1;
}

int GEN_BANDES_12::NED3_Build_Check_Morphs() {
	// clear band3 redundancy
	if (wb1b2.ntab1b2 && wb1b2.Check3(&grid0[54], &grid0[54])) return 0;	
	tww.MorphToP1(*this);
	n_auto_p1 = tblnauto[it16_3]; //ia 0-415 not index
	t_auto_p1 = &automorphsp[tblnautostart[it16_3]];
	if (n_auto_p1) {// build min canonical Bx
		wb1b2.Init(grid1);
		int zs1[81]; memcpy(zs1, grid1, sizeof zs1);
		// move to mini
		for (int imorph = 0; imorph < n_auto_p1; imorph++) {
			BANDMINLEX::PERM& p = t_auto_p1[imorph];
			int* z = &grid1[27];// morph the band
			SKT_MORPHTOP
				BandReOrder(band);
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
		memcpy(grid1, zs1, sizeof grid1);
	}
	if (irtw2	// check direct against min
		&&tww2.IsBelowGrid1(*this, tww1.zs0, 2))return 0;	
	if (i2t16 == i1t16) {// grid1 Bx2=Bx3
		if (grid1[27] != 1) return 0;
		if (tww.IsBelowGrid1(*this, grid1, 0))return 0; // acb
	}
	if (i3t16 == i1t16 	// grid1 Bx1=Bx2
		&& tww.IsBelowGrid1(*this, grid1, 1))return 0;   // bac
	return 1;
}
int GEN_BANDES_12::F3B_See_NED3() {//  pass 1 c<=a<=b 

	if (!NED3_Build_Check_Morphs()) return 0;
	if (i3t16 < i1t16 && i1t16 < i2t16) return 2; 
	if (i2t16 == i3t16) {// same 3 bands  
		for (int i = 0; i <= 5; i++)
			if (tww.IsBelowGrid1(*this, grid1, i))return 0;
		if (irtw2)for (int i = 0; i < 6; i++) {
			if (tww.IsBelowGrid1(*this, tww1.zs0, i))return 0;
		}
		return 4;
	}
	if (irtw2 == 2) {// see diagonal
		if (i1t16 == i3t16 &&  // do 
			tww2.IsBelowGrid1(*this, tww1.zs0, 0))return 0;
		if (i1t16 == i2t16 && // do cab perm cba
			tww2.IsBelowGrid1(*this, tww1.zs0, 3))return 0;
	}
	if (i1t16 == i3t16) 	return 8;	
	return 1;
}


 