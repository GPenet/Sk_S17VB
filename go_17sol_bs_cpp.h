#ifdef XXXXXX
//___ start process expand bands collect uas guas ...
const char * diagband =  "364912578715834962892675341"; // la bonne
const char * diagband3 = "276341895538297416941568237";
const char * diagpuz = "...45..........6...89...1....4......7.5.........6..3.......1.95.3.2.............7";
#endif

#define TEST_ON
//#define CHECKHARVEST



G17B::G17B() {
	gang27 = gang[0];
	debug17 = debug17_check=aigstop = 0;
	// Init SGUAs tables load permanent data
	int tpcol[3][2] = { {1,2},{0,2},{0,1} };
	for (int ist = 0; ist < 3; ist++) {//stack
		for (int ircol = 0; ircol < 3; ircol++) {//first column
			int* p = tpcol[ircol], dcol = 3 * ist;
			for (int irdig = 0; irdig < 3; irdig++) {//digit1
				for (int irdig2 = 0; irdig2 < 3; irdig2++) {//digit 2
					int i = 27 * ist + 9 * ircol + 3 * irdig + irdig2;
					SGUA2& w = tsgua2[i];
					w.i_81 = i; w.i9= dcol + ircol;
					w.col1 = dcol + p[0];		w.col2 = dcol + p[1];
					w.id1 = irdig + 3 * w.col1;
					w.id2 = irdig2 + 3 * w.col2;
					if (0) {
						cout << "sua2=" << w.i_81 //<< " i9=" << w.i9 + 1
							<< " cols12 " << w.col1 + 1 << w.col2 + 1
							<< " id1;id2 " << w.id1 << ";" << w.id2 << endl;
					}
				}
			}
		}
	}
	for (int ist = 0; ist < 3; ist++) {//stack
		for (int irdig1 = 0; irdig1 < 3; irdig1++) {//digit1
			for (int irdig2 = 0; irdig2 < 3; irdig2++) {//digit 2
				for (int irdig3 = 0; irdig3 < 3; irdig3++) {//digit 2
					int i = 27 * ist + 9 * irdig1 + 3 * irdig2 + irdig3;
					SGUA3& w = tsgua3[i];
					w.i_81 = i;
					w.stack = ist;
					w.col1 = 3 * ist;// minirow first column in gangster
					// pointers to gangster digits
					w.id1 = 9 * ist+irdig1;
					w.id2 = 9 * ist + 3+irdig2;
					w.id3 = 9 * ist + 6 + irdig3;
					if (0) {
						cout << "sgua3 " << w.i_81 << " col1=" << w.col1 + 1
							<< " id1;id2,id3 " << w.id1
							<< ";" << w.id2 << ";" << w.id3 << endl;
					}
				}
			}
		}
	}
}
void G17B::StartInit() {
	memcpy(grid0, genb12.grid0, sizeof grid0);// use first b3
	memcpy(&grid0[54], genb12.bands3[0].band0, sizeof genb12.bands3[0].band0);
	zh2b_g.Init_g0(grid0);// init brute force bands 1+2
	if (sgo.vx[4]==2)mincluesb3 = 7;
	else mincluesb3 = 6;
	// Build gangster 3x3
	for (int i = 0; i < 9; i++) {// 9 cols to build out of gangcols
		int istack = C_box[i];
		int* d = gang[i], c = genb12.gangcols[i];
		uint32_t bit;
		bitscanforward(bit, c);
		c ^= (1 << bit);
		d[0] = bit;
		gang_digits_cols[bit][istack] = i;
		bitscanforward(bit, c);
		c ^= (1 << bit);
		d[1] = bit;
		gang_digits_cols[bit][istack] = i;
		bitscanforward(bit, c);
		d[2] = bit;
		gang_digits_cols[bit][istack] = i;
	}

	//finish sguas2/3 set up with the gangster 
	gsock2.SetAll_0();  gsock3.SetAll_0();
	//register int* gang27=genb12.gang27; // redefines gang[9][3] as 27 integer
	//for (int i = 0; i < 27; i++) cout << gang27[i] + 1;
	//cout << endl;
	int tpcol[3][2] = { {1,2},{0,2},{0,1} };
	for (register int ist = 0,i=0,stack=07007007; ist < 3; ist++,stack<<=3) {
		//cout << Char27out(stack) << endl;
		for (register int irst = 0; irst < 27; irst++, i++) {//stack
			{//________________________ guas2
				SGUA2& w = tsgua2[i];
				w.dig1 = gang27[w.id1];
				w.dig2 = gang27[w.id2];
				w.digs = (1 << w.dig1) | (1 << w.dig2);
				w.valid = 0;
				//cout << i << " " << w.dig1+1 <<" "<< w.dig2+1 << endl;
				//continue;
				for (register int ib3 = 0; ib3 < genb12.nband3; ib3++) {
					STD_B3& b3 = genb12.bands3[ib3];
					register uint32_t Bf = b3.dpat[w.dig1] | b3.dpat[w.dig2],
						Bfc = (Bf | (Bf >> 9) | (Bf >> 18)) & 0777;// colummn pattern
					Bf &= stack;
					int nr = 0, bfr = 0, irow;
					if (Bf & 0777) { nr++; bfr |= Bf & 0777; irow = 1; }// row1
					if (Bf & 0777000) { nr++; bfr |= Bf & 0777000; irow = 2; }
					if (Bf & 0777000000) { nr++; bfr |= Bf & 0777000000; irow = 3; }
					if (nr == 1) {// we have a gua2
						//cout << ib3 << " " << i << endl;
						gsock2.setBit(i);
						w.valid = 1;
						b3.g.gsocket2.setBit(i);
						b3.g.pat2[i] = Bf;
						int imini = 3 * irow + ist, mask = 7 << (3 * imini);
						b3.g.ua2_imini[i] = imini;
						int bit27= mask ^ Bf,i27;
						bitscanforward(i27, bit27);
						b3.g.ua2pair27[i] = bit27;
						b3.i_27_to_81[i27] = i;
					}
					else {// gua46 except if 6 columns
						int ncol = _popcnt32(Bfc);
						if (ncol == 5) {// gua2_4
							b3.g.pat2[i] = bfr;
						}
						else b3.g.pat2[i] = Bf;
					}
				}
			}
			
			
			{ //________________________ guas3

				SGUA3& w = tsgua3[i];
				w.dig1 = gang27[w.id1];
				w.dig2 = gang27[w.id2];
				w.dig3 = gang27[w.id3];
				w.digs = 1 << w.dig1 | 1 << w.dig2 | 1 << w.dig3;
				w.valid = 0;
				for (register int ib3 = 0; ib3 < genb12.nband3; ib3++) {
					STD_B3& b3 = genb12.bands3[ib3];
					register uint32_t Bf = b3.dpat[w.dig1]
						| b3.dpat[w.dig2] | b3.dpat[w.dig3];// colummn pattern
					Bf &= stack;
					int nr = 0,  irow=0;
					if (Bf & 0777) { nr++;  irow =0; }// row1
					if (Bf & 0777000) { nr++;  irow = 1; }
					if (Bf & 0777000000) { nr++;  irow = 2; }
					if (nr == 1) {// we have a gua3
						gsock3.setBit(i);
						w.valid = 1;
						b3.g.gsocket3.setBit(i);
						b3.g.pat3[i] = Bf;
						int imini = 3 * irow + ist;
						b3.g.ua3_imini[i] = imini;
						b3.i_9_to_81[imini] = i;
					}
				}
			}			
		}
	}
}
void G17B::StartInitDebug() {
	cout << "check guas2 status" << endl;
	for (register int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		STD_B3& b3 = genb12.bands3[ib3];
		cout << Char64out(b3.g.gsocket2.bf.u64[0]);
		cout << Char27out(b3.g.gsocket2.bf.u32[2])<<" "
			<< b3.g.gsocket2 .Count96() << endl;
	}
	cout << Char64out(gsock2.bf.u64[0]);
	cout << Char27out(gsock2.bf.u32[2]) << " all guas2 "
		<< gsock2.Count96() << endl;

	cout << "check guas3 status" << endl;
	for (register int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		STD_B3& b3 = genb12.bands3[ib3];
		cout << Char64out(b3.g.gsocket3.bf.u64[0]);
		cout << Char27out(b3.g.gsocket3.bf.u32[2]) 
			<< b3.g.gsocket3.Count96() << endl;
	}
	cout << Char64out(gsock3.bf.u64[0]);
	cout << Char27out(gsock3.bf.u32[2]) << " all guas3 " 
		<< gsock3.Count96() << endl;

}
void G17B::StartPrint() {
	is_test_on = 1;
	cout << myband2.band << " go band2 id=" << myband2.i416 << " nb12=" << genb12.nb12
		<< " nb3=" << genb12.nband3 << " p_cpt2g[0]=" << p_cpt2g[0] << endl;
	cout << Char64out(gsock2.bf.u64[0]);
	cout << Char27out(gsock2.bf.u32[2]) << " all guas2 "
		<< gsock2.Count96() << endl;
	cout << Char64out(gsock3.bf.u64[0]);
	cout << Char27out(gsock3.bf.u32[2]) << " all guas3 "
		<< gsock3.Count96() << endl;
	if (0) {
		cout << myband2.band << endl << endl;
		for (int i = 0; i < genb12.nband3; i++) {
			STD_B3& b3 = genb12.bands3[i];
			cout << b3.band << " b3 i=" << i << endl;
			//for (int id = 0; id < 9; id++)
			//	cout << Char27out(b3.dpat[id]) << endl;
		}
	}
	for (int i = 0; i < 81; i++) cout << grid0[i] + 1;
//	cout << " grid to use in ua collector" << endl;

	//StartInitDebug();

}

void G17B::Start() {// processing an entry 
	if (aigstop)return;
#ifdef TEST_ONEB2
	if (p_cpt2g[0]) return;
#endif
	is_test_on = 0;
	p_cpt2g[0] ++;
	p_cpt2g[1] += genb12.nband3;

	// back if min 66 or not 5 clues in band 2  pass 3

	StartInit();// do all preliminary setups
#ifdef TEST_ON
	StartPrint();
#endif
	UaCollector();//_do uas/guas2_3 initial harvest
	StartAfterUasHarvest();
}
#ifdef DEBUGKNOWN

void G17B::GoM10Known() {// processing an entry 656 566 with the relevant table of ba,ds3
//if(1)return;
	if (aigstop)return;
	p_cpt2g[0] ++;
	if (g17b.debug17 > 1) diag = diagbug = g17b.debug17;
	else diag = diagbug = 0;

}
#else
void G17B::StartKnown() {
}
#endif
 
void G17B::UaCollector() {
	FirstUasCollect();// uas 2 digits bands and stacks
	//return;
	SecondUasCollect();// uas 345 bands 1+2
	UasCollect6_7();
	Guas2Collect();
	UasCollect4box();// 4 box in bands 1+2
}

inline void G17B::Adduab12( uint32_t digs, uint32_t nd) {

	uint64_t* t = zh2gxn.tua;
	uint32_t n = zh2gxn.nua,n2=0;
	BF128 tw[100];
	for (uint32_t i = 0; i < n; i++) {
		register uint64_t w = t[i], cc = _popcnt64(w);
		if (cc < 25) {	
			BF128 x; x.bf.u64[0] =  w; x.bf.u64[1] = cc;	
			tw[n2++] = x;
		}
	}
	if (n2 > 1) {// sort to check redundancy
		for (uint32_t i = 0; i < n2 - 1; i++) {
			for (uint32_t j=i+1; j < n2 ; j++) {
				if (tw[i].bf.u64[1] > tw[j].bf.u64[1]) {
					BF128 x = tw[i]; tw[i] = tw[j]; tw[j] = x;
				}
			}
		}
	}
	for (uint32_t i = 0; i < n2; i++) {
		BF128 x = tw[i];
		register uint64_t cc = x.bf.u64[1],
			w= x.bf.u64[0],nw=~w;

		if (i) {
			for (uint32_t j = 0; j < i; j++) {
				BF128 y = tw[j];
				register uint64_t ycc = y.bf.u64[1];
				if (ycc == cc) break;
				if(! (y.bf.u64[0] & nw)) { cc = 0; break; }
			}
			if (!cc)continue;// subset found
		}
		if (cc > 15) cc = 15;
		w |= (cc << 59);
		if (tuasb12.nua < UA12SIZE)tuasb12.AddInit(w, digs, nd);
	}

}

void G17B::FirstUasCollect() {// produce all uas 2/3 digits

	struct TUAS81 {// used to collect uas 2 digits
		BF128 	tall[200];// first collect
		uint32_t ntall;// , ntold;
		int Add(BF128 w, uint32_t floor) {
			uint32_t cc = w.Count96(), nfloor = ~floor;
			for (uint32_t iua = 0; iua < ntall; iua++) {
				BF128 wt = tall[iua];
				uint32_t cct = wt.bf.u32[3] >> 16;
				uint32_t floort = wt.bf.u32[3] & 0777;
				wt.bf.u32[3] = 0;
				if (cct < cc) {// look for subset
					if (floort & nfloor) continue;//can not be subset
					if ((wt - w).isEmpty()) return 0; // subset
					continue;
				}
				if (wt == w) return 0; // redundancy
				if (cct == cc)continue; // not sorted here
				// insert here and check super set
				for (uint32_t jua = ntall; jua > iua; jua--)
					tall[jua] = tall[jua - 1];
				tall[iua] = w;// new inserted
				tall[iua].bf.u32[3] = floor | (cc << 16);
				ntall++;
				// is it a subset of a previous entry
				for (iua++; iua < ntall; iua++) {
					if ((w - tall[iua]).isEmpty()) {// we have a subset
						for (uint32_t k = iua + 1; k < ntall; k++)
							tall[k - 1] = tall[k];
						ntall--;
						iua--; //continue same position
					}
				}
				return 2;
			}
			w.bf.u32[3] = floor | (cc << 16);
			tall[ntall++] = w;// new added
			return 1;
		}

		void DebugAll() {
			cout << "all uas status nuas=" << ntall << endl;
			for (uint32_t iua = 0; iua < ntall; iua++) {
				BF128 wt = tall[iua];
				cout << Char27out(wt.bf.u32[0]) << "|";
				cout << Char27out(wt.bf.u32[1]) << "|";
				cout << Char27out(wt.bf.u32[2]) << " size=";
				cout << wt.bf.u16[7]
					<< " i=" << iua
					<< " " << _popcnt64(wt.bf.u64[0])
					<< " " << _popcnt32(wt.bf.u32[2]) << endl;
			}
		}

	}tuas81;
	tuas81.ntall = 0;// no ua so far
	zhgxn.SetupFsol(grid0);
	for (int i = 0; i < 36; i++) {
		uint32_t myf = floors_2d[i];
		zhou2[0].GoZ2(myf);
		for (uint32_t j = 0; j < zhgxn.nua; j++) {
			BF128 w = zhgxn.tua[j];
			tuas81.Add(w, myf);
		}
	}
	/*
	// insert bands  and apply subsets
	uint32_t ntsort[27 - 4];
	memset(ntsort, 0, sizeof ntsort); // empty table
	{ //sort all by size
		BF128* t = tuas81.tall;
		uint32_t nt = tuas81.ntall;
		for (uint32_t iua = 0; iua < nt; iua++) {
			BF128 wt = t[iua];
			uint32_t cc = wt.bf.u16[7] - 4;// index is 0 for 4
			wt.bf.u32[3] = 0; // clear digits and count
			w_tfua[cc][ntsort[cc]++] = wt;
		}
	}
	// here insert all missing uas one band (1 or 2)
	for (int i = 0; i < 2; i++) {
		STD_B1_2 & wb=(i)? myband2: myband1;
		BF128 wt;
		wt.SetAll_0();
		uint32_t* tu = wb.tua, nu = wb.nua;
		for (uint32_t j = 0; j < nu; j++) {
			register uint32_t u = tu[j] & BIT_SET_27, cc = _popcnt32(u);
			cc -= 4;
			wt.bf.u32[i] = u;
			w_tfua[cc][ntsort[cc]++] = wt;
		}
	}
	// reload tall and check subsets/redundancy
	tuas81.ntall = 0;
	for (int i = 0; i < 23; i++)if (ntsort[i]) {
		BF128* tt = w_tfua[i];
		for (uint32_t j = 0; j < ntsort[i]; j++) {
			tuas81.Add2(tt[j], i + 4);
		}
	}
	*/
	// split uas band 12 and others
	chunkh.nt2digs	= tuasb12.nua = 0;
	for (uint32_t iua = 0; iua < tuas81.ntall; iua++) {
		BF128 wt = tuas81.tall[iua];
		if (wt.bf.u32[2]) 	chunkh.t2digs[chunkh.nt2digs++] = wt;
		else {
			uint64_t w = wt.bf.u64[0], cc = _popcnt64(w);
			if(cc > 15 )cc = 15;
			//cout << Char64out(w) << " " << cc << " w cc 64 + ??" << endl;
			w |= cc << 59;
			tuasb12.AddInit(w, wt.bf.u16[6],2);
		}
	}
	//tuas81.DebugAll();
	//tuasb12.DumpInit();
}

struct EXTUAS {// get a sorted subset of uas for next "collect" step
	uint64_t  t2[300];
	uint32_t  nt2;
	void GetSortB(uint64_t* t, uint32_t n,uint64_t filter, uint32_t nold = 0) {
		BF128 vsize[25][4];
		uint64_t  t2a[512];
		uint32_t  nt2a=0;
		nt2 = nold;
		register uint64_t Fn = (~filter) & BIT_SET_2X, w;
		memset(vsize, 0, sizeof vsize);
		for (uint32_t iua = 0; iua < tuasb12.nua; iua++) {
			w = tuasb12.tua[iua];
			uint64_t cc = w >> 59;
			w &= BIT_SET_2X;
			if (!(w & Fn)) {
				int bloc = nt2a >> 7, ir = nt2a - (bloc << 7);
				vsize[cc][bloc].setBit(ir);
				t2a[nt2a++] = w;// clean the count 
			}
			if (nt2a >= 512)break;
		}
		//cout << nt2a << " nt2a" << endl;
		uint32_t nbl64 = (nt2a + 63) >> 6, x;
		for (int i1 = 4; i1 < 25; i1++) {
			uint64_t* tb64 = vsize[i1]->bf.u64;
			for (uint32_t i2 = 0; i2 < nbl64; i2++) if (tb64[i2]) {
				register uint64_t V = tb64[i2];
				while (bitscanforward64(x, V)) {// switch to 54 mode here
					V ^= (uint64_t)1 << x;
					t2[nt2++] = t2a[x + (i2 << 6)];
					if (nt2 >= 100)return;

				}
			}
		}		


	}
	void GetSortA(uint64_t filter, uint32_t nold = 0) {
		nt2 = nold;
		register uint64_t Fn = (~filter) & BIT_SET_2X, w;
		for (uint32_t iua = 0; iua < tuasb12.nua; iua++) {
			w = tuasb12.tua[iua];
			if (!(w & Fn))t2[nt2++] = w;
			if (nt2 >= 100)break;
		}
		if (nt2) {
			if (nt2 > 1)
				for (uint32_t i = 0; i < nt2 - 1; i++) {
					register uint64_t Ri = t2[i];
					for (uint32_t j = i + 1; j < nt2; j++) {
						register uint64_t  Rj = t2[j];
						if ((Ri >> 59) > (Rj >> 59)) { // smaller size first 
							t2[i] = Rj; t2[j] = Ri; Ri = Rj;
						}
					}
				}
			for (uint32_t i = 0; i < nt2; i++)
				t2[i] &= BIT_SET_2X;// clean the count

		}
	}
	void Dump() {
		cout << "dump extuas nt2=" << nt2 << endl;
		for (uint32_t i = 0; i < nt2; i++)
			cout << Char2Xout(t2[i]) << "i=" << i << endl;
	}
}extuas;

void G17B::SecondUasCollect() {// collect 345 digits in bands 1+2
	//cout << "seconuacollect" << endl;
	//zh2gxn.InitKnown(tuasb12.t2, &tuasb12.nt2);
	zh2gxn.InitKnown(extuas.t2, &extuas.nt2);
	zh2gxn.SetupFsol(grid0);
	//cout << "3 digits" << endl;
	for (int i = 0; i < 84; i++) {
		uint32_t myf = floors_3d[i];
		int ir = zh2_3[0].GoZ3(myf);// cells unsolved count
		if (ir < 6) continue;// minimum for a fresh ua 3 digits
		uint64_t F = zh2_3[0].cells_unsolved.bf.u64;
		extuas.GetSortA(F);
		zh2_3[0].DoZ3Go();
		if (zh2gxn.nua) Adduab12(myf,3);
	}
	if (is_test_on)tuasb12.DumpShort("3 digits");
	for (int i = 0; i < 126; i++) {
		int ir = zh2_4[0].GoZ4(floors_4d[i]);
		if (ir < 8) continue;// minimum for a fresh ua 4 digits
		uint64_t F = zh2gxn.unsolved_field;
		extuas.GetSortA(F);
		zh2_4[0].DoZ4Go();
		if (zh2gxn.nua) Adduab12(floors_4d[i],4);
	}
	if (is_test_on)tuasb12.DumpShort("4 digits");
	for (int i = 0; i < 126; i++) {
		uint32_t myf = 0777 ^ floors_4d[i];
		int ir = zh2_5[0].GoZ5(myf);
		uint64_t F = zh2gxn.unsolved_field;
		//extuas.GetSort(F);
		extuas.GetSortB(tuasb12.tua, tuasb12.nua, F);
		zh2_5[0].DoZ5Go();
		if (zh2gxn.nua) Adduab12(myf, 5);
	}
	if (is_test_on)tuasb12.DumpShort("5 digits");
}
void G17B::UasCollect6_7() {
	for (int i = 0; i < 84; i++) {
		uint32_t ass= floors_3d[i], myf = 0777^ass;
		uint64_t bf = zh2gxn.Getsol(ass);
		uint64_t F = bf^BIT_SET_2X;
		//cout << Char2Xout(bf) << " ";cout << Char9out(myf) << endl;
		extuas.GetSortB(tuasb12.tua, tuasb12.nua, F);
		//extuas.GetSort(F);
		zh2b[0].InitBf(bf);
		zh2b[0].Dob();
		if (zh2gxn.nua) Adduab12(myf, 6);
	}
	if (is_test_on)tuasb12.DumpShort("6 digits");
	return;// 7 seems too expensive 
	for (int i = 0; i < 36; i++) {
		uint32_t ass = floors_2d[i], myf = 0777 ^ ass;
		uint64_t bf = zh2gxn.Getsol(ass);
		uint64_t F = bf ^ BIT_SET_2X;
		extuas.GetSortB(tuasb12.tua, tuasb12.nua,F);
		zh2b[0].InitBf(bf);
		zh2b[0].Dob();
		if (zh2gxn.nua) Adduab12(myf, 7);
	}
}
void G17B::UasCollect4box() {
	//cout << "4 box collect" << endl;
	//___________________ box 1245
	{
		uint64_t stack = units3xBM[5].u64[0], b4 = BIT_SET_2X & (~stack);
		extuas.GetSortB(tuasb12.tua, tuasb12.nua, b4);
		//extuas.Dump();
		zh2b[0].InitBf(stack);

	}
	zh2b[0].Dob(22);
	if (zh2gxn.nua) {
		for (uint32_t i = 0; i < zh2gxn.nua; i++) {
			register uint64_t w = zh2gxn.tua[i], cc = _popcnt64(w);	
			//cout << Char2Xout(w) << " to add b4 size " << cc << endl;
			if (tuasb12.nua < 2550)tuasb12.AddInit(w | (cc << 59), 0,0);
		}
	}
	
	//___________________ box 1346
	{
		uint64_t stack = units3xBM[4].u64[0], b4 = BIT_SET_2X & (~stack);
		extuas.GetSortB(tuasb12.tua, tuasb12.nua,b4);
		//extuas.Dump();
		zh2b[0].InitBf(stack);

	}
	zh2b[0].Dob(22);
	if (zh2gxn.nua) {
		for (uint32_t i = 0; i < zh2gxn.nua; i++) {
			register uint64_t w = zh2gxn.tua[i], cc = _popcnt64(w);
			//cout << Char2Xout(w) << " to add b4 size " << cc << endl;
			if (tuasb12.nua < UA12SIZE)tuasb12.AddInit(w | (cc << 59), 0,0);
		}
	}
	//___________________ box 2356
		{
		uint64_t stack = units3xBM[5].u64[3], b4 = BIT_SET_2X & (~stack);
		//extuas.GetSort(b4);
		extuas.GetSortB(tuasb12.tua, tuasb12.nua, b4);
		//extuas.Dump();
		zh2b[0].InitBf(stack);

	}
	zh2b[0].Dob(22);
	if (zh2gxn.nua) {
		for (uint32_t i = 0; i < zh2gxn.nua; i++) {
			register uint64_t w = zh2gxn.tua[i], cc = _popcnt64(w);	
			//cout << Char2Xout(w) << " to add b4 size " << cc << endl;
			if (tuasb12.nua < UA12SIZE)tuasb12.AddInit(w | (cc << 59), 0,0);
		}
	}
	if(is_test_on)tuasb12.DumpShort("4 box ");
}

struct TTT1 {//tuab12;tua <=12 per size, then per floor
#define NB 3
	uint64_t tua0[NB*128], tua[NB * 128];
	uint32_t nua, tdigs0[NB * 128], tdigs[NB * 128],nbl;
	BF128 vs [9][NB] ;// vectors size 4/12
	struct VF {
		BF128 v2[36][NB], v3[84][NB], v4[126][NB], v5[126][NB], v6[84][NB];
	}vf;
	void Init(uint64_t* tu, uint32_t* td, uint32_t n) {
		//cout << "uas <=12 out of n="<<n << endl;
		memset(vs, 0, sizeof vs);
		nua = 0;
		for (uint32_t i = 0; i < n; i++) {
			register uint64_t U = tu[i],cc=U>>59;
			if (cc>=4 && cc <= 12) {
				tua0[nua] = U & BIT_SET_2X;
				tdigs0[nua] = td[i];
				uint32_t bloc = nua / 128, ir = nua - 128 * bloc;
				vs[cc-4][bloc].setBit(ir);
				nua++;
				if (nua >= NB * 128) break;// overflow control
			}
		}
		//cout << "uas <=12 got nua=" << nua << endl;
		// copy per size 
		nbl = (nua + 63) / 64;
		uint32_t x, n2 = 0;
		for (int i = 0; i < 9; i++) {// size 4/12
			uint64_t* v = vs[i][0].bf.u64;
			for (uint32_t j = 0; j < nbl; j++) {// uas for this size
				uint64_t V = v[j];
				while (bitscanforward64(x, V)) {
					V ^= (uint64_t)1 << x;
					register uint32_t ii = x + 64 * j;
					tua[n2] = tua0[ii]  ;
					tdigs[n2++] = tdigs0[ii];
				}
			}
		}
		//DumpUas();
		// build vectors per floor
		memset(&vf, 0, sizeof vf);
		for (int i = 0; i < 36; i++) {// ______________2 digits
			register uint32_t fl = floors_2d[i];
			BF128* pv = vf.v2[i];
			for (uint32_t j = 0; j < nua; j++)
				if (tdigs[j] == fl) {
					uint32_t bloc =j/ 128, ir = j - 128 * bloc;
					pv[bloc].setBit(ir);
				}
		}
		for (int i = 0; i < 84; i++) {// ______________3 digits
			register uint32_t fl = floors_3d[i],fln=~fl;
			BF128* pv = vf.v3[i];
			for (uint32_t j = 0; j < nua; j++)
				if (tdigs[j] == fl) {
					uint32_t bloc = j / 128, ir = j - 128 * bloc;
					pv[bloc].setBit(ir);
				}
			for (int i = 0; i < 36; i++) {// subset___2 digits
				uint32_t fl2 = floors_2d[i];
				if (fl2 & fln) continue;
				BF128* pv2 = vf.v2[i];
				for (uint32_t k = 0; k < nbl; k++) pv[k] |= pv2[k];
			}
		}
		for (int i = 0; i < 126; i++) {// ______________4 digits
			register uint32_t fl = floors_4d[i], fln = ~fl;
			BF128* pv = vf.v4[i];
			for (uint32_t j = 0; j < nua; j++)
				if (tdigs[j] == fl) {
					uint32_t bloc = j / 128, ir = j - 128 * bloc;
					pv[bloc].setBit(ir);
				}
			for (int i = 0; i < 84; i++) {// subset___3 digits
				uint32_t fl2 = floors_3d[i];
				if (fl2 & fln) continue;
				BF128* pv2 = vf.v3[i];
				for (uint32_t k = 0; k < nbl; k++) pv[k] |= pv2[k];
			}
		}
		for (int i = 0; i < 126; i++) {// ______________5 digits
			register uint32_t fl = 0777^floors_4d[i], fln = ~fl;
			BF128* pv = vf.v5[i];
			for (uint32_t j = 0; j < nua; j++)
				if (tdigs[j] == fl) {
					uint32_t bloc = j / 128, ir = j - 128 * bloc;
					pv[bloc].setBit(ir);
				}
			for (int i = 0; i < 126; i++) {// subset___4 digits
				uint32_t fl2 = floors_4d[i];
				if (fl2 & fln) continue;
				BF128* pv2 = vf.v4[i];
				for (uint32_t k = 0; k < nbl; k++) pv[k] |= pv2[k];
			}
		}
		for (int i = 0; i < 84; i++) {// ______________6 digits
			register uint32_t fl = 0777 ^ floors_3d[i], fln = ~fl;
			BF128* pv = vf.v6[i];
			for (uint32_t j = 0; j < nua; j++)
				if (tdigs[j] == fl) {
					uint32_t bloc = j / 128, ir = j - 128 * bloc;
					pv[bloc].setBit(ir);
				}
			for (int i = 0; i < 126; i++) {// subset___5 digits
				uint32_t fl2 = 0777^floors_4d[i];
				if (fl2 & fln) continue;
				BF128* pv2 = vf.v5[i];
				for (uint32_t k = 0; k < nbl; k++) pv[k] |= pv2[k];
			}
		}
		//DumpV();
	}
	int Load(uint64_t* td, int iv, int ix) {
		BF128* v;
		switch (iv) {
		case 2: v = vf.v2[ix]; break;
		case 3: v = vf.v3[ix]; break;
		case 4: v = vf.v4[ix]; break;
		case 5: v = vf.v5[ix]; break;
		case 6: v = vf.v6[ix]; break;
		default: return 0;
		}
		int n = 0,x;
		uint64_t* pv = v->bf.u64;
		for (uint32_t i = 0; i < nbl; i++) {
			uint64_t V = pv[i];
			while (bitscanforward64(x, V)) {
				V ^= (uint64_t)1 << x;
				register uint32_t ii = x + 64 * i;
				td[n++] = tua[ii];
			}
		}
		return n;
	}
	void DumpUas() {
		cout << "uas <=12 sorted" << endl;
		if (nua < 256) for (uint32_t i = 0; i < nua; i++) {
			cout << Char2Xout(tua[i]) << "|";
			cout  << Char9out(tdigs[i]) << " i="<<i << endl;
		}
	}
	void DumpV() {
		cout << " dump v nbl=" << nbl << endl;
		cout << "v2" << endl;
		for (int i = 0; i < 36; i++) {
			uint64_t *pv=vf.v2[i]->bf.u64;
			if (!(*pv)) continue;
			for (uint32_t k = 0; k < nbl; k++)
				cout << Char64out(pv[k])<< "|";
			cout << Char9out(floors_2d[i]) <<" i="<<i << endl;

		}
		cout << "v3" << endl;
		for (int i = 0; i < 84; i++) {
			uint64_t* pv = vf.v3[i]->bf.u64;
			if (!(*pv)) continue;
			for (uint32_t k = 0; k < nbl; k++)
				cout << Char64out(pv[k])<< "|";
			cout << Char9out(floors_3d[i]) << " i=" << i << endl;

		}
		cout << "v4" << endl;
		for (int i = 0; i < 126; i++) {
			uint64_t* pv = vf.v4[i]->bf.u64;
			if (!(*pv)) continue;
			for (uint32_t k = 0; k < nbl; k++)
				cout << Char64out(pv[k]) << "|";
			cout << Char9out(floors_4d[i]) << " i=" << i << endl;

		}
		cout << "v5" << endl;
		for (int i = 0; i < 126; i++) {
			uint64_t* pv = vf.v5[i]->bf.u64;
			if (!(*pv)) continue;
			for (uint32_t k = 0; k < nbl; k++)
				cout << Char64out(pv[k]) << "|";
			cout << Char9out(0777^floors_4d[i]) << " i=" << i << endl;

		}
		cout << "v6" << endl;
		for (int i = 0; i < 84; i++) {
			uint64_t* pv = vf.v6[i]->bf.u64;
			if (!(*pv)) continue;
			for (uint32_t k = 0; k < nbl; k++)
				cout << Char64out(pv[k]) << "|";
			cout << Char9out(0777^floors_3d[i]) << " i=" << i << endl;

		}

	}
	//36;84;126;126 for 2;3;4;5; 

}ttt1;

struct GUAH {// handler for initial collection of guas 2 3

	struct GUA {
		uint64_t tua[128];
		uint32_t nua, type, i81;
		inline void Add(uint64_t ua) {
			if (nua < 128)tua[nua++] = ua;
		}
		inline void Init(uint32_t n, uint32_t t, uint32_t i) {
			nua = n; type = t; i81 = i; 
		}
		int Load(uint64_t* tu, uint64_t bf) {
			int n=0;
			register uint64_t nbf = ~bf;
			for (uint32_t i = 0; i < nua; i++) {
				register uint64_t U = tua[i],cc=_popcnt64(U);
				if (n> 10) break;// 
				if (n > 5 &&cc>12) break;// 
				if (!(U & nbf))tu[n++] = U;
			}
			return n;
		}
		void SortClean() {
			if (nua < 2)return;
			GUA w = *this;
			BF128 vsize[30];
			memset(vsize, 0, sizeof vsize);
			uint64_t tcopy[128];
			memcpy(tcopy, w.tua, sizeof tcopy);
			for (uint32_t i = 0; i < nua; i++) {
				register uint64_t  cc = _popcnt64(tua[i]&BIT_SET_2X);
				vsize[cc].setBit(i);
			}
			nua = 0;
			for (int i = 0; i < 20; i++) if (vsize[i].isNotEmpty()) {
				int j;
				BF128 v = vsize[i];
				while ((j = v.getFirst128()) >= 0) {
					v.clearBit(j);
					register uint64_t U = tcopy[j];
					for (uint32_t k = 0; k < nua; k++) {
						register uint64_t U2 = tua[k];
						if (!(U & (~U2))) { U = 0; break; }
					}
					if (U) tua[nua++] = U;
					if (nua > 40) break;
				}
				if (nua > 40) break;
			}

		}
		void Debug(int nodet = 1) {
			cout << "gua type=" << type << " i81=" << i81
				<< "\tnua=" << nua << endl;
			if (nodet) return;
			for (uint32_t i = 0; i < nua; i++) {
				cout << Char2Xout(tua[i])<< " " << _popcnt64(tua[i] & BIT_SET_2X)
					<< " " << i << endl;
			}
		}
	}tg2[81],tg3[81],guaw;
	void Init() {
		for (int i = 0; i < 81; i++) {
			tg2[i].Init(0, 0, i); tg3[i].Init(0, 1, i);
		}
	}
	inline void Add2(uint64_t u, uint32_t i) { 	if(_popcnt64(u)<18)		tg2[i].Add(u); }
	inline void Add3(uint64_t u, uint32_t i) { if (_popcnt64(u) < 18)		tg3[i].Add(u); }

	int IsUa4(int i81) {
		GUA& g = tg2[i81];
		if (g.nua == 1 && _popcnt64(g.tua[0]) == 2) return 1;
		return 0;
	}
	int IsUamin(int i81) {
		GUA& g = tg3[i81];
		if (g.nua == 1 && _popcnt64(g.tua[0]) < 5) return 1;
		return 0;
	}

	void SortClean3() {
		for (int i = 0; i < 81; i++)
			if (tg3[i].nua > 1) tg3[i].SortClean();
	}
	void SortClean() {
		for (int i = 0; i < 81; i++)
			if (tg2[i].nua > 1) tg2[i].SortClean();
	}

	int CutG2(int lim) {
		int n = 0;
		for (int i = 0; i < 81; i++) {
			if ((int)tg2[i].nua > lim)tg2[i].nua = lim;
			n += tg2[i].nua;
		}
		return n;
	}
	int CutG3(int lim) {
		int n = 0;
		for (int i = 0; i < 81; i++) {
			if ((int)tg3[i].nua > lim)tg3[i].nua = lim;
			n += tg3[i].nua;
		}
		return n;
	}
	void Dump2all2() {
		for (int i = 0; i < 81; i++) if(tg2[i].nua) {
			tg2[i].Debug(0);
		}
	}
	void Dump2all3() {
		for (int i = 0; i < 81; i++) if (tg3[i].nua) {
			tg3[i].Debug(0);
		}
	}
}guah;
struct GUAH54 {// handler guas 2 3 in 54 mode
	uint64_t gbuf[162 * 40]; // room for cut 30+10 in average
	struct GUA {
		uint64_t *tua, killer;
		uint32_t nua, nuamax, type, i81;
		inline void Init(uint64_t *p, uint32_t t, uint32_t i) {
			tua=p; type = t; i81 = i; killer = ~0;
			nua = nuamax = 0;
		}

		void Debug(int nodet = 1) {
			cout << "gua type=" << type << " i81=" << i81
				<< "\tnua=" << nua << endl;
			cout  << Char54out(killer) << " K" << endl;
			if (nodet) return;
			for (uint32_t i = 0; i < nua; i++) {
				cout << Char54out(tua[i]) << " " << _popcnt64(tua[i] )
					<< " " << i << endl;
			}
		}
	}tg2[81], tg3[81];
	void Build();

	void Dumpall2() {
		for (int i = 0; i < 81; i++)  {
			tg2[i].Debug(0);
		}
	}
	void Dumpall3() {
		for (int i = 0; i < 81; i++)  {
			tg3[i].Debug(0);
		}
	}
}guah54;
void GUAH54::Build() {// cut to 30 switch to 54 killer

	uint64_t* pbuf = gbuf;
	for (int i81 = 0; i81 < 81; i81++){
		tg2[i81].Init(pbuf, 0, i81);
		if (g17b.gsock2.On(i81)) {
			GUAH::GUA& g0 = guah.tg2[i81];
			GUA& gd = tg2[i81];
			uint32_t n = g0.nua;
			if (n > 30) n = 30;
			gd.nua = n;
			gd.nuamax = gd.nua + 10;
			pbuf += gd.nuamax;// lock the storing place
			register uint64_t K = ~0;
			for (uint32_t i = 0; i < n; i++) {
				register uint64_t U = g0.tua[i];
				U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54
				K &= U;
				gd.tua[i] = U;
			}
			gd.killer = K;
		}
	}
	for (int i81 = 0; i81 < 81; i81++) {
		tg3[i81].Init(pbuf, 1, i81);
		if (g17b.gsock3.On(i81)) {
			GUAH::GUA& g0 = guah.tg3[i81];
			GUA& gd = tg3[i81];
			uint32_t n = g0.nua;
			if (n > 30) n = 30;
			gd.nua = n;
			gd.nuamax = gd.nua + 10;
			pbuf += gd.nuamax;// lock the storing place
			register uint64_t K = ~0;
			for (uint32_t i = 0; i < n; i++) {
				register uint64_t U = g0.tua[i];
				U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54
				K &= U;
				gd.tua[i] = U;
			}
			gd.killer = K;
		}
	}
	//Dumpall2();
	//Dumpall3();
}


int HasUaHitFalse(uint64_t* tu, uint32_t nu0, uint32_t nu, uint32_t i81) {
	SGUA2 w = tsgua2[i81];
	register uint64_t bf = zh2b_g.fd_sols[0][w.dig1].bf.u64;
	bf |= zh2b_g.fd_sols[0][w.dig2].bf.u64;
	uint64_t bf2 = units3xBM[9 + w.col1].u64[0];
	bf2 |= units3xBM[9 + w.col2].u64[0];
	bf &= bf2;
	for(uint32_t i=nu0;i<nu;i++) if(bf&tu[i]) return 1;
	return 0;
}
void G17B::Guas2Collect() {
	ttt1.Init(tuasb12.tua, tuasb12.tdigs, tuasb12.nua);
	guah.Init();

	for (int i81 = 0; i81 < 81; i81++) if (gsock2.On(i81)) {// initial gua2s 4 cells
		SGUA2 w = tsgua2[i81];
		uint64_t bf = zh2b_g.fd_sols[0][w.dig1].bf.u64;
		bf |= zh2b_g.fd_sols[0][w.dig2].bf.u64;
		uint64_t bf2 = units3xBM[9 + w.col1].u64[0];
		bf2 |= units3xBM[9 + w.col2].u64[0];
		bf &= bf2;
		int xcell1, xcell2;
		bitscanforward64(xcell1, bf);
		bitscanreverse64(xcell2, bf);
		int cell1=From_128_To_81[xcell1], cell2= From_128_To_81[xcell2];
		int r1 = C_row[cell1], r2 = C_row[cell2];
		if (r1 == r2) {
			//cout << Char54out(bf) << " gua2_4 i_81=" << i81 << endl;
			guah.Add2(bf, i81);
		}
	}
	/*
	chunkh.Init();
	cur_ib = 0;
	for (uint32_t iua = 0; iua < chunkh.nt2digs; iua++) {
		BF128 w = chunkh.t2digs[iua];
		if (w.bf.u64[0]) {
			int cc = BuildGua(w);
			if (cc == 2) {
				cout << Char54out(w.bf.u64[0]) << " i_81=" << w.bf.u32[2] << endl;
				guah.Add2(w.bf.u64[0], w.bf.u32[2]);
			}
			//chunkh.Add128(w, cc);
		}
		//else chunkh.Addband3(w.bf.u32[2]);//redundancy
	}
	*/




	Guas2CollectG2();
	Guas2CollectG3();

}
void G17B::Guas2CollectG2() {
	//___ find guas2 2 digits
	for (int i81 = 0; i81 < 81; i81++) if (gsock2.On(i81)) {
		if (guah.IsUa4(i81)) continue;// nothing to do
		SGUA2 w = tsgua2[i81];
		GUAH::GUA& gt = guah.tg2[i81];
		int i = 0;
		for (i; i < 36; i++)if (floors_2d[i] == w.digs) break;
		extuas.nt2 = ttt1.Load(extuas.t2, 2, i);
		if (extuas.nt2 && HasUaHitFalse(extuas.t2,0, extuas.nt2, i81))
			continue;

		int ir = zh2_2[0].GoZ2G2(w.digs, w.col1, w.dig1, w.col2, w.dig2);
		if (ir < 0)continue; // locked
		if (ir == 1) {//	solved)
			uint64_t U = zh2gxn.tua[0], cc = _popcnt64(U);
			guah.Add2(U, i81);
		}
		uint64_t F = zh2gxn.unsolved_field;
		zh2_2[0].DoZ2Go();
		if (zh2gxn.nua)
			for (uint32_t i = 0; i < zh2gxn.nua; i++) {
				uint64_t U = zh2gxn.tua[i];
				guah.Add2(U, i81);
			}
	}
	guah.SortClean();
	//___ find guas2 3 digits

	for (int i81 = 0; i81 < 81; i81++) if (gsock2.On(i81)) {
		if (guah.IsUa4(i81)) continue;// nothing to do
		SGUA2 w = tsgua2[i81];
		GUAH::GUA& gt = guah.tg2[i81];
		for (int i = 0; i < 84; i++) {// find UAs 3 digits
			int fl = floors_3d[i];
			if (!((fl & w.digs) == w.digs)) continue;
			uint64_t isfl = 0;// cells ok of the floors
			for (int i2 = 0, bit = 1; i2 < 9; i2++, bit <<= 1)
				if (fl & bit)	isfl |= zh2gxn.fsol[i2];
			extuas.nt2=gt.Load(extuas.t2, isfl);
			zh2gxn.nkguas = extuas.nt2;
			uint32_t istart = extuas.nt2;
			extuas.nt2 += ttt1.Load(&extuas.t2[istart], 3, i);
			// check no ua hitting cells "false"
			if (extuas.nt2>istart && HasUaHitFalse(extuas.t2, istart,extuas.nt2, i81))
				continue;
			int ir = zh2_3[0].GoZ3G2(fl, w.col1, w.dig1, w.col2, w.dig2),ir2;
			if (ir < 0)continue; // locked
			if (ir == 1) {//	solved)
				if (zh2gxn.nua) {
					uint64_t U = zh2gxn.tua[0];
					guah.Add2(U, i81);
				}
				continue;
			}
			ir2 = zh2_3[0].DoZ3Go();
			if (ir2 < 0) continue;
			if (zh2gxn.nua)
				for (uint32_t i = 0; i < zh2gxn.nua; i++) {
					uint64_t U = zh2gxn.tua[i];
					guah.Add2(U, i81);
				}
		}
	}// end 3 digits
	//___ find guas2 4 digits
	guah.SortClean();

	for (int i81 = 0; i81 < 81; i81++) if (gsock2.On(i81)) {
		if (guah.IsUa4(i81)) continue;// nothing to do
		SGUA2 w = tsgua2[i81];
		GUAH::GUA& gt = guah.tg2[i81];
		for (int i = 0; i < 126; i++) { 
			int fl = floors_4d[i];
			if (!((fl & w.digs) == w.digs)) continue;
			uint64_t isfl = 0;// cells ok of the floors
			for (int i2 = 0, bit = 1; i2 < 9; i2++, bit <<= 1)
				if (fl & bit)	isfl |= zh2gxn.fsol[i2];
			extuas.nt2 = gt.Load(extuas.t2, isfl);
			zh2gxn.nkguas = extuas.nt2;
			uint32_t istart = extuas.nt2;
			extuas.nt2 += ttt1.Load(&extuas.t2[istart], 3, i);
			// check no ua hitting cells "false"
			if (extuas.nt2 > istart && HasUaHitFalse(extuas.t2, istart, extuas.nt2, i81))
				continue;
			int ir = zh2_4[0].GoZ4G2(fl, w.col1, w.dig1, w.col2, w.dig2), ir2;
			if (ir < 0)continue; // locked
			if (ir == 1) {//	solved)
				if (zh2gxn.nua) {
					uint64_t U = zh2gxn.tua[0];
					guah.Add2(U, i81);
				}
				continue;
			}
			ir2 = zh2_4[0].DoZ4Go();
			if (ir2 < 0) continue;
			if (zh2gxn.nua)
				for (uint32_t i = 0; i < zh2gxn.nua; i++) {
					uint64_t U = zh2gxn.tua[i];
					guah.Add2(U, i81);
				}
		}
	}// end 4 digits

	guah.SortClean();

	//___ find guas2 5 digits

	for (int i81 = 0; i81 < 81; i81++) if (gsock2.On(i81)) {
		if (guah.IsUa4(i81)) continue;// nothing to do
		SGUA2 w = tsgua2[i81];
		GUAH::GUA& gt = guah.tg2[i81];
		for (int i = 0; i < 126; i++) { 
			int fl =0777^ floors_4d[i];
			if (!((fl & w.digs) == w.digs)) continue;
			uint64_t isfl = 0;// cells ok of the floors
			for (int i2 = 0, bit = 1; i2 < 9; i2++, bit <<= 1)
				if (fl & bit)	isfl |= zh2gxn.fsol[i2];
			extuas.nt2 = gt.Load(extuas.t2, isfl);
			zh2gxn.nkguas = extuas.nt2;
			uint32_t istart = extuas.nt2;
			extuas.nt2 += ttt1.Load(&extuas.t2[istart], 3, i);
			// check no ua hitting cells "false"
			if (extuas.nt2 > istart && HasUaHitFalse(extuas.t2, istart, extuas.nt2, i81))
				continue;
			int ir = zh2_5[0].GoZ5G2(fl, w.col1, w.dig1, w.col2, w.dig2), ir2;
			if (ir < 0)continue; // locked
			if (ir == 1) {//	solved)
				if (zh2gxn.nua) {
					uint64_t U = zh2gxn.tua[0];
					guah.Add2(U, i81);
				}
				continue;
			}
			ir2 = zh2_5[0].DoZ5Go();
			if (ir2 < 0) continue;
			if (zh2gxn.nua)
				for (uint32_t i = 0; i < zh2gxn.nua; i++) {
					uint64_t U = zh2gxn.tua[i];
					guah.Add2(U, i81);
				}
		}
	}// end 5 digits
	guah.SortClean();
	//___ find guas2 6 digits
	for (int i81 = 0; i81 < 81; i81++) if (gsock2.On(i81)) {
		if (guah.IsUa4(i81)) continue;// nothing to do
		SGUA2 w = tsgua2[i81];
		GUAH::GUA& gt = guah.tg2[i81];
		if (gt.nua >= 30) continue;
		for (int i = 0; i < 84; i++) {
			int fl = 0777 ^ floors_3d[i];
			if (!((fl & w.digs) == w.digs)) continue;
	
			uint64_t isfl = 0;// cells ok of the floors
			for (int i2 = 0, bit = 1; i2 < 9; i2++, bit <<= 1)
				if (fl & bit)	isfl |= zh2gxn.fsol[i2];
			extuas.nt2 = gt.Load(extuas.t2, isfl);
			zh2gxn.nkguas = extuas.nt2;
			uint32_t istart = extuas.nt2;
			extuas.nt2 += ttt1.Load(&extuas.t2[istart], 3, i);
			// check no ua hitting cells "false"
			if (extuas.nt2 > istart && HasUaHitFalse(extuas.t2, istart, extuas.nt2, i81))
				continue;
			zh2b[0].InitBfG2(BIT_SET_2X^ isfl, w.col1, w.dig1, w.col2, w.dig2);
			//continue;
			zh2b[0].Dob(14);
			if (zh2gxn.nua)
				for (uint32_t i = 0; i < zh2gxn.nua; i++) {
					uint64_t U = zh2gxn.tua[i];
					//if (_popcnt64(U) >= 16) continue;
					//cout <<Char2Xout(	U) << " added6 i81=" << i81 
					//	<< "  "<< _popcnt64(U) << endl;
					guah.Add2(U, i81);
				}
		}
	}// end 6 digits

	guah.SortClean();
	ng2=guah.CutG2(30);
	if (is_test_on)cout << "ng2 after cut30=" << ng2 << endl;
	//guah.Dump2all();

}
void G17B::Guas2CollectG3() {
	for (int i81 = 0; i81 < 81; i81++) if (gsock3.On(i81)) {
		//if (i81)continue;
		SGUA3 w = tsgua3[i81];
		GUAH::GUA& gt = guah.tg3[i81];
		// Setup the perms for gangsters in minirow
		int bita = 1 << w.dig1, bitb = 1 << w.dig2, bitc = 1 << w.dig3,
			digs = w.digs,triplet_perms[2][3];

		int* p = triplet_perms[0];// for per abc -> bca
		p[0] = bita | bitb; p[1] = bitb | bitc; p[2] = bitc | bita;

		p = triplet_perms[1];// for per abc -> cab
		p[0] = bita | bitc; p[1] = bitb | bita; p[2] = bitc | bitb;
		int tp3f[2][3] = { {1,2,0},{2,0,1} };// perms no valid digit
		for (int ip = 0; ip < 2; ip++) {
			// build revised gangster
			int rgangcols[9];// revised gangster
			memcpy(rgangcols, zh2gxn.gangsters, sizeof rgangcols);
			p = triplet_perms[ip];
			int c1 = w.col1, c2 = c1 + 1, c3 = c1 + 2;
			rgangcols[c1] ^= p[0];
			rgangcols[c2] ^= p[1];
			rgangcols[c3] ^= p[2];
			int i = 0;
			for (i; i < 84; i++)if (floors_3d[i] == w.digs) break;
			extuas.nt2 = ttt1.Load(extuas.t2, 3, i);
			//extuas.Dump();
			// find UAs 3 digits one 3 digits 
			int ir= zh2_3[0].GoZ3G3(w.digs, rgangcols),ir2;
			if (ir < 0)continue; // locked
			//if (extuas.nt2 && HasUaHitFalse3(extuas.t2, 0, extuas.nt2, i81))
			//	continue;
			//cout << "i81=" << i81 << " digs=" << w.dig1 + 1 << w.dig2 + 1 << w.dig3 + 1
			//	<< " col1=" << c1 + 1<<" ir ="<<ir << endl;
			//zh2_3[0].ImageCandidats();
			if (ir == 1) {//	solved)
				if (zh2gxn.nua) {
					uint64_t U = zh2gxn.tua[0];
					//cout << Char2Xout(U) << "to add " <<_popcnt64(U) << endl;
					guah.Add3(U, i81);
				}
				continue;
			}
			ir2 = zh2_3[0].DoZ3Go();
			if (ir2 < 0) continue;
			if (zh2gxn.nua)
				for (uint32_t i = 0; i < zh2gxn.nua; i++) {
					uint64_t U = zh2gxn.tua[i];
					//cout << Char2Xout(U) << "to add " << _popcnt64(U) << endl;
					guah.Add3(U, i81);
				}
		}
	}
	//guah.Dump2all3();
	Guas2CollectG3_4d();
}
void G17B::Guas2CollectG3_4d() {
	for (int i81 = 0; i81 < 81; i81++) if (gsock3.On(i81)) {
		SGUA3 w = tsgua3[i81];
		GUAH::GUA& gt = guah.tg3[i81];
		if (guah.IsUamin(i81)) continue;// nothing to do

		// Setup the perms for gangsters in minirow
		int bita = 1 << w.dig1, bitb = 1 << w.dig2, bitc = 1 << w.dig3,
			digs = w.digs, triplet_perms[2][3];

		int* p = triplet_perms[0];// for per abc -> bca
		p[0] = bita | bitb; p[1] = bitb | bitc; p[2] = bitc | bita;

		p = triplet_perms[1];// for per abc -> cab
		p[0] = bita | bitc; p[1] = bitb | bita; p[2] = bitc | bitb;
		int tp3f[2][3] = { {1,2,0},{2,0,1} };// perms no valid digit
		for (int ip = 0; ip < 2; ip++) {
			//cout << i81 << " perm " << ip << endl;
			// build revised gangster
			int rgangcols[9];// revised gangster
			memcpy(rgangcols, zh2gxn.gangsters, sizeof rgangcols);
			p = triplet_perms[ip];
			int c1 = w.col1, c2 = c1 + 1, c3 = c1 + 2;
			rgangcols[c1] ^= p[0];
			rgangcols[c2] ^= p[1];
			rgangcols[c3] ^= p[2];
			for (int i = 0; i < 126; i++) {
				int fl = floors_4d[i];
				if (!((fl & w.digs) == w.digs)) continue;
				uint64_t isfl = 0;// cells ok of the floors
				for (int i2 = 0, bit = 1; i2 < 9; i2++, bit <<= 1)
					if (fl & bit)	isfl |= zh2gxn.fsol[i2];
				extuas.nt2 = gt.Load(extuas.t2, isfl);
				zh2gxn.nkguas = extuas.nt2;
				uint32_t istart = extuas.nt2;
				extuas.nt2 += ttt1.Load(&extuas.t2[istart], 4, i);
				//cout << "nolds=" << istart << " ";
				int ir = zh2_4[0].GoZ4G3(fl, rgangcols), ir2;
				if (ir < 0)continue; // locked
				//extuas.Dump();
				//zh2_4[0].ImageCandidats();
				if (ir == 1) {//	solved)
					if (zh2gxn.nua) {
						uint64_t U = zh2gxn.tua[0];
						if (_popcnt64(U) < 18) {
							//cout << Char2Xout(U) << "to add " << _popcnt64(U) << endl;
							guah.Add3(U, i81);
						}
					}
					continue;
				}
				ir2 = zh2_4[0].DoZ4Go();
				if (ir2 < 0) continue;
				if (zh2gxn.nua)
					for (uint32_t i = 0; i < zh2gxn.nua; i++) {
						uint64_t U = zh2gxn.tua[i];
						if (_popcnt64(U) < 18) {
							//cout << Char2Xout(U) << "to add " << _popcnt64(U) << endl;
							guah.Add3(U, i81);
						}
					}
			}
		}
	}
	guah.SortClean3();
	//guah.Dump2all3();
	Guas2CollectG3_5d();
}
void G17B::Guas2CollectG3_5d() {
	for (int i81 = 0; i81 < 81; i81++) if (gsock3.On(i81)) {
		SGUA3 w = tsgua3[i81];
		GUAH::GUA& gt = guah.tg3[i81];
		if (guah.IsUamin(i81)) continue;// nothing to do
		if (gt.nua >= 10) continue;

		// Setup the perms for gangsters in minirow
		int bita = 1 << w.dig1, bitb = 1 << w.dig2, bitc = 1 << w.dig3,
			digs = w.digs, triplet_perms[2][3];

		int* p = triplet_perms[0];// for per abc -> bca
		p[0] = bita | bitb; p[1] = bitb | bitc; p[2] = bitc | bita;

		p = triplet_perms[1];// for per abc -> cab
		p[0] = bita | bitc; p[1] = bitb | bita; p[2] = bitc | bitb;
		int tp3f[2][3] = { {1,2,0},{2,0,1} };// perms no valid digit
		for (int ip = 0; ip < 2; ip++) {
			//cout << i81 << " perm " << ip << endl;
			// build revised gangster
			int rgangcols[9];// revised gangster
			memcpy(rgangcols, zh2gxn.gangsters, sizeof rgangcols);
			p = triplet_perms[ip];
			int c1 = w.col1, c2 = c1 + 1, c3 = c1 + 2;
			rgangcols[c1] ^= p[0];
			rgangcols[c2] ^= p[1];
			rgangcols[c3] ^= p[2];
			for (int i = 0; i < 126; i++) {
				int fl = 0777^floors_4d[i];
				if (!((fl & w.digs) == w.digs)) continue;
				uint64_t isfl = 0;// cells ok of the floors
				for (int i2 = 0, bit = 1; i2 < 9; i2++, bit <<= 1)
					if (fl & bit)	isfl |= zh2gxn.fsol[i2];
				extuas.nt2 = gt.Load(extuas.t2, isfl);
				zh2gxn.nkguas = extuas.nt2;
				uint32_t istart = extuas.nt2;
				extuas.nt2 += ttt1.Load(&extuas.t2[istart], 5, i);
				int ir = zh2_5[0].GoZ5G3(fl, rgangcols), ir2;
				//cout << "nolds=" << istart << " nt2=" << extuas.nt2 << "ir=" <<ir << endl;

				if (ir < 0)continue; // locked
				//extuas.Dump();
				//zh2_4[0].ImageCandidats();
				if (ir == 1) {//	solved)
					if (zh2gxn.nua) {
						uint64_t U = zh2gxn.tua[0];
						if (_popcnt64(U) < 18) {
							//cout << Char2Xout(U) << "to add " << _popcnt64(U) << endl;
							guah.Add3(U, i81);
						}
					}
					continue;
				}

				ir2 = zh2_5[0].DoZ5Go();
				if (ir2 < 0) continue;
				if (zh2gxn.nua)
					for (uint32_t i = 0; i < zh2gxn.nua; i++) {
						uint64_t U = zh2gxn.tua[i];
						if (_popcnt64(U) < 18) {
							//cout << Char2Xout(U) << "to add " << _popcnt64(U) << endl;
							guah.Add3(U, i81);
						}
					}
			}
		}
	}
	guah.SortClean3();
	ng3 = guah.CutG3(20);
	if (is_test_on)cout << "ng3 after cut20=" << ng3 << endl;
	//guah.Dump2all3();
}

//__________ end of initial uas harvest switch to 54 mode

void G17B::StartAfterUasHarvest() {
#ifdef CHECKHARVEST
	if (p_cpt2g[0] == 1) {// min,max,count,b1 max
		p_cpt2g[33] = p_cpt2g[38] = p_cpt2g[43] = p_cpt2g[48] = 0;
		p_cpt2g[30] = p_cpt2g[31] = p_cpt2g[32] = genb12.nband3;
		p_cpt2g[35] = p_cpt2g[36] = p_cpt2g[37] = tuasb12.nua;
		p_cpt2g[40] = p_cpt2g[41] = p_cpt2g[42] = ng2;
		p_cpt2g[45] = p_cpt2g[46] = p_cpt2g[47] = ng3;
		return;
	}
	int nx= genb12.nband3;
	if (nx < p_cpt2g[30])p_cpt2g[30] = nx;
	if (nx > p_cpt2g[31]) {p_cpt2g[31] = nx; p_cpt2g[33] = p_cpt2g[0];}
	p_cpt2g[32] += nx;

	nx = tuasb12.nua;
	if (nx < p_cpt2g[35])p_cpt2g[35] = nx;
	if (nx > p_cpt2g[36]) { p_cpt2g[36] = nx; p_cpt2g[38] = p_cpt2g[0]; }
	p_cpt2g[37] += nx;

	if (ng2 < p_cpt2g[40])p_cpt2g[40] = ng2;
	if (ng2 > p_cpt2g[41]) { p_cpt2g[41] = ng2; p_cpt2g[43] = p_cpt2g[0]; }
	p_cpt2g[42] += ng2;

	if (ng3 < p_cpt2g[45])p_cpt2g[45] = ng3;
	if (ng3 > p_cpt2g[46]) { p_cpt2g[46] = ng3; p_cpt2g[48] = p_cpt2g[0]; }
	p_cpt2g[47] += ng3;
	return;
#endif
	guah54.Build();
	//guah.Setupvectors();
	//guah.CheckSetup();
	//return;
	t54b12.Build_ta128(tuasb12.tua, tuasb12.nua);
	//t54g2.Build_ta128();
	//t54g3.Build_ta128();
	Expand_03();
	//tuasb12.DumpInit();

	//PutUasStartInVector();

}

void T54B12::Build_ta128(uint64_t* t, uint32_t n) {
	InitA();
	BF128 vsize[25][30];
	uint32_t nbl64 = (n + 63) >> 6, x;
	memset(vsize, 0, sizeof vsize);
	for (uint32_t i = 0; i < n; i++) {
		int bloc = i >> 7, ir = i - (bloc << 7);
		uint64_t cc = t[i] >> 59;
		vsize[cc][bloc].setBit(ir);
	}

	uint64_t tw[UA12SIZE]; // to check redundancy
	uint32_t ntw = 0;
	for (int i1 = 4; i1 < 25; i1++) {
		uint64_t* tb64 = vsize[i1]->bf.u64;
		for (uint32_t i2 = 0; i2 < nbl64; i2++) if (tb64[i2]) {
			register uint64_t V = tb64[i2];
			while (bitscanforward64(x, V)) {// switch to 54 mode here
				V ^= (uint64_t)1 << x;
				register uint64_t R = t[x + (i2 << 6)];
				R = (R & BIT_SET_27) | ((R & BIT_SET_B2) >> 5);
				if (1) {// check redundancy and subsets
					for (uint32_t i = 0; i < ntw; i++) {
						if (!(R & (~tw[i]))) {
							cout << Char54out(R) << " erased" << endl;;
							cout << Char54out(tw[i]) <<" i="<<i << endl;
							R = 0; break;						}
					}
					if (R)tw[ntw++] = R;
					else continue;
				}
				AddA(R);
			}
		}
	}
	cout << "ntw=" << ntw << " nold" << n << endl;
	if(ntw!=n)ta128[0].Dump();
}


//________ start expand uas bands 12

void G17B::Expand_03() {
	SPB03* s, * sn;
	if (aigstop) return;
	//zh2b[0].InitBands12(grid0);
	T54B12::TUVECT& tuv128 = t54b12.ta128[0];
	uint64_t* twu = tuv128.t;
	s = spb_0_15;	memset(s, 0, sizeof spb_0_15[0]);
	s->active_cells = maskLSB[54].u64[0];
	s->possible_cells = twu[0];
	s->v = tuv128.v0;// initial nothing done

	//_______ start search 3 first clues
next:
	// catch and apply cell in bitfields
	int cell;
	uint64_t p = s->possible_cells;
	if (!p)	if (--s >= spb_0_15)goto next; else return;
	bitscanforward64(cell, p);
	register uint64_t bit = (uint64_t)1 << cell;
	s->possible_cells ^= bit;
	tclues[s->ncl] = cell;
	s->active_cells ^= bit;
	uint64_t ac = s->active_cells;
	sn = s + 1; *sn = *s; sn->ncl++;
	sn->cbs.Add(cell);
	sn->all_previous_cells |= bit;
	sn->v &= tuv128.vc[cell];
	if (sn->ncl == 3) {// 3 cellsfirst step
		p_cpt2g[3]++;
#ifdef HAVEKNOWN
		if (!((~pk54) & sn->all_previous_cells)) {
			cout << Char54out(sn->all_previous_cells) << " expected 3" << endl;
			Expand_46();
			aigstop = 1;
		}
		else 	Expand_46();
#else
		if (_popcnt64(ac) < 45) {
			cout << Char54out(sn->all_previous_cells) << " 3clues [3] " << p_cpt2g[3] << endl;
			if (t54b12.Build_tb128()) goto next;
			Expand_46();
		}
#endif
		if (aigstop) return;
		goto next;
	}
	// find next ua
	int ir = sn->v.getFirst128();
	if (ir < 0) return;//never
	uint64_t Ru = twu[ir] & ac;
	if (!Ru)goto next;//dead branch unlikely
	sn->possible_cells = Ru;
	s++; // switch to next spot
	goto next;
}


int T54B12::Build_tb128() {
	BF128 tvw[20];
	uint32_t lastbloc = t54b12.nablocs;
	tvw[0] = spb_0_15[3].v;
	for (uint32_t i = 1; i <= lastbloc; i++) {
		TUVECT & vv= ta128[i];
		BF128 v = vv.v0 , * vc=vv.vc;
		for (uint32_t ic = 0; ic < 3; ic++)
			v &= vc[g17b.tclues[ic]];
		tvw[i]=v;
	}

	// apply active on still valid uas and flag by size
	BF128 vsize[19][UABNBLOCS];
	memset(vsize, 0, sizeof vsize);
	//uint32_t nbl64 = (n + 63) >> 6, x;
	uint64_t tw[128* UABNBLOCS]; 
	uint32_t ntw = 0;
	{
		register uint64_t Ac = spb_0_15[3].active_cells;
		for (uint32_t i = 0; i <= lastbloc; i++) {
			register uint64_t * t= ta128[i].t;
			BF128 V= tvw[i];
			while (1) {
				register int ir = V.getFirst128();
				if (ir>=0) {
					V.clearBit(ir);
					register uint64_t R = t[ir] & Ac;
					if (!R)return 1; //dead
					register uint64_t cc = _popcnt64(R);
					if (cc > 18)cc = 18;
					int bloc = ntw >> 7, ir = ntw - (bloc << 7);
					vsize[cc][bloc].setBit(ir);
					tw[ntw++] = R;
					if (ntw >= 128 * UABNBLOCS) break;
				}
				else break;
			}
			if (ntw >= 128 * UABNBLOCS) break;
		}
	}
	//__Build the reduced set of UAs vectors
	InitB();
	uint32_t nbl64 = (ntw + 63) >> 6, x;
	for (int i1 = 1; i1 < 19; i1++) {
		uint64_t* tb64 = vsize[i1]->bf.u64;
		for (uint32_t i2 = 0; i2 < nbl64; i2++) if (tb64[i2]) {
			register uint64_t V = tb64[i2];
			while (bitscanforward64(x, V)) {
				V ^= (uint64_t)1 << x;
				AddB(tw[x + (i2 << 6)]);
			}
		}
	}
	if (g17b.is_test_on)cout << "ua3c ntw=" << ntw << " nb128=" << nb128 << endl;

	return 0;
}

void G17B::Expand_46() {
	if (aigstop) return;
	SPB03*sl= &spb_0_15[4] ,* s=sl, * sn;
	T54B12::TUVECT& tu128 = t54b12.tb128[0];
	uint64_t* twu = tu128.t;
	*s = spb_0_15[3];	// duplicate 3 for new vector
	s->possible_cells = twu[0];
	s->v = tu128.v0;// initial nothing done

	//_______ start search clues 4-6
next:	// catch and apply cell in bitfields
	register int cell;
	uint64_t p = s->possible_cells;
	if (!p)	if (--s >= sl)goto next; else return;
	bitscanforward64(cell, p);
	register uint64_t bit = (uint64_t)1 << cell;
	s->possible_cells ^= bit;
	tclues[s->ncl] = cell;
	s->active_cells ^= bit;
	uint64_t ac = s->active_cells;
	sn = s + 1; *sn = *s; sn->ncl++;
	sn->cbs.Add(cell);
	sn->all_previous_cells |= bit;
	sn->v &= tu128.vc[cell];
	if (sn->ncl == 6) {// 6 cells
		p_cpt2g[4]++;
		if (_popcnt64(ac) < 36) {
			cout << Char54out(sn->all_previous_cells) << "\t 6clues [4] " << p_cpt2g[4] 
				<<" "<< _popcnt64(ac) << endl;
			if (t54b12.Build_tc128()) goto next;
			t54g2.BuildG2();
			t54g3.BuildG3();
		}

		//SetupExpand_7p();
#ifdef HAVEKNOWN
		if (!((~pk54) & sn->all_previous_cells)) {
			cout << Char54out(sn->all_previous_cells) << " expected 6" << endl;
			if (mincluesb3 == 6)Expand_7_11();
			else 		Expand_7_10();
			aigstop = 1;
		}
		else {
			if (mincluesb3 == 6)Expand_7_11();
			else 		Expand_7_10();
		}
#else
		//if (mincluesb3 == 6)Expand_7_11();
		//else 		Expand_7_10();
#endif
		if (aigstop) return;
		goto next;
	}
	// find next ua
	int ir = sn->v.getFirst128();
	uint64_t Ru;
	if (ir < 0) {//never valid
		if (t54b12.nbblocs) {// more uas to check
			for (uint32_t i = 1; i <= t54b12.nbblocs; i++) {
				T54B12::TUVECT& vv = t54b12.tb128[i];
				BF128 v = vv.v0, * vc = vv.vc;
				for (uint32_t ic = 3; ic < sn->ncl; ic++)
					v &= vc[tclues[ic]];
				if (v.isNotEmpty()) {
					int ir2 = v.getFirst128();
					uint64_t Ru = vv.t[ir2] & ac;
					if (!Ru)goto next;//dead branch
					sn->possible_cells = Ru;
					s++; // switch to next spot
					goto next;
				}
			}

		}
		if (zh2b[1].IsValid(tclues, sn->ncl)) {
			uint32_t i = zh2gxn.nua - 1;
			register uint64_t ua = zh2gxn.tua[i],
				cc = _popcnt64(ua),
				ua54 = (ua & BIT_SET_27) | ((ua & BIT_SET_B2) >> 5);
			t54b12.AddA(ua54);
			t54b12.AddB(ua54);
			Ru = ua;
		}
		else {
			cout << "bug exp 4-7 lower 7" << endl;
			aigstop = 1; return;
		}

	}
	else  Ru = twu[ir] & ac;
	if (!Ru)goto next;//dead branch
	sn->possible_cells = Ru;
	s++; // switch to next spot
	goto next;
}



int T54B12::Build_tc128() {
	BF128 tvw[20];
	uint32_t lastbloc = t54b12.nbblocs;
	tvw[0] = spb_0_15[7].v;
	for (uint32_t i = 1; i <= lastbloc; i++) {
		TUVECT& vv = tb128[i];
		BF128 v = vv.v0, * vc = vv.vc;
		for (uint32_t ic = 3; ic < 6; ic++)
			v &= vc[g17b.tclues[ic]];
		tvw[i] = v;
	}

	// apply active on still valid uas and flag by size
	BF128 vsize[19][UACNBLOCS];
	memset(vsize, 0, sizeof vsize);
	uint64_t tw[128 * UACNBLOCS]; 
	uint32_t ntw = 0;
	{
		register uint64_t Ac = spb_0_15[7].active_cells;
		for (uint32_t i = 0; i <= lastbloc; i++) {
			register uint64_t* t = tb128[i].t;
			BF128 V = tvw[i];
			while (1) {
				register int ir = V.getFirst128();
				if (ir >= 0) {
					V.clearBit(ir);
					register uint64_t R = t[ir] & Ac;
					if (!R)return 1; //dead
					register uint64_t cc = _popcnt64(R);
					if (cc > 18)cc = 18;
					int bloc = ntw >> 7, ir = ntw - (bloc << 7);
					vsize[cc][bloc].setBit(ir);
					tw[ntw++] = R;
					if (ntw >= 128 * UACNBLOCS) break;
				}
				else break;
			}
			if (ntw >= 128 * UACNBLOCS) break;
		}
	}
	//__Build the reduced set of UAs vectors clean redundancy
	InitC();
	uint32_t nbl64 = (ntw + 63) >> 6, x;
	for (int i1 = 1; i1 < 19; i1++) {
		uint64_t* tb64 = vsize[i1]->bf.u64;
		for (uint32_t i2 = 0; i2 < nbl64; i2++) if (tb64[i2]) {
			register uint64_t V = tb64[i2];
			while (bitscanforward64(x, V)) {
				V ^= (uint64_t)1 << x;
				register uint64_t U = tw[x + (i2 << 6)];
				// check redundancy in tc128[0]
				if(t54b12.IsNotRedundant(U))
					AddC(U);
			}
		}
	}
	if (g17b.is_test_on)cout << "ua6c ntw=" << ntw << " nc128=" << nc128 << endl;
	return 0;
}
void T54G2::BuildG2() {//extract, shrink, sort g2
	Init();
	BF128 vsize[25][30];
	memset(vsize, 0, sizeof vsize);
	uint64_t tw[128 * 30];// can accept 40 uas per i81
	uint32_t twi[128 * 30], twn[128 * 30]; // i81,size 
	uint32_t ntw = 0;
	{ // extract still valid
		register uint64_t F = spb_0_15[7].all_previous_cells,
			A = spb_0_15[7].active_cells;
		cout << Char54out(F) << "F " << endl;
		cout << Char54out(A) << "A " << endl;
		for (uint32_t i81 = 0; i81 < 81; i81++) {
			GUAH54::GUA& gt = guah54.tg2[i81];
			if (gt.killer & F)continue; // killed or not valid
			uint32_t n1 = gt.nua;
			//cout << "i81=" << i81 << " n1=" << n1 << endl;
			if (n1 == 1) {// killer is enough
				register uint64_t U = gt.tua[0]&A;
				twi[ntw] = i81;
				twn[ntw] = (uint32_t)_popcnt64(U);
				tw[ntw++] = U;
				continue;
			}

			// must go in details and clean redundant
			if (n1 > 40) n1 = 40;
			uint64_t tw2[40];// valid
			uint32_t ntw2 = 0,iw2;
			uint64_t vsize2[25];// sorting tw2 by size
			memset(vsize2, 0, sizeof vsize2);

			// store still valid for active cells
			for (uint32_t i = 0; i < n1; i++) {
				register uint64_t U = gt.tua[i];
				if (U & F)continue; // hit
				U &= A;// keep active cells
				uint64_t cc = _popcnt64(U);
				vsize2[cc] |= (uint64_t)1 << ntw2;
				tw2[ntw2++] = U;
			}
			if (vsize2[0]) {// force true (one not hit)
				twi[ntw] = i81;
				twn[ntw] = 0;
				tw[ntw++] = 0;
				continue;
			}
			// take back per size and send to main table
			uint32_t istart = ntw; // start for redundancy
			for (int i = 0; i < 20; i++)  {
				register uint64_t V = vsize2[i];
				while (bitscanforward64(iw2, V)) {
					V ^= (uint64_t)1 << iw2;
					register uint64_t U = tw2[iw2],	nU=~U;
					for(uint32_t j=istart;j<ntw;j++)
						if (!(tw[j] & nU)) { U = 0; break; }
					if (U) {
						twi[ntw] = i81;
						twn[ntw] = i;
						tw[ntw++] = U;
					}
				}
			}
		}// end i81 
	}// end build tw
	//cout << " buildg2 ntw=" << ntw << endl;


	// sort tw by size 
	for (uint32_t i = 0; i < ntw; i++) {
		uint32_t bloc = i >> 7,ir = i - (bloc << 7);
		vsize[twn[i]][bloc].setBit(ir);
	}
	// send the sorted table tw to t54g2
	cout << "BuildG2 back g2 ntw=" << ntw << endl;

	uint32_t nbl64 = (ntw + 63) >> 6,x,n=0;

	for (int i1 = 0; i1 < 25; i1++) {
		uint64_t* tb64 = vsize[i1]->bf.u64;
		for (uint32_t i2 = 0; i2 < nbl64; i2++) if (tb64[i2]) {
			register uint64_t V = tb64[i2];
			while (bitscanforward64(x, V)) {
				V ^= (uint64_t)1 << x;
				uint32_t xx = x + (i2 << 6);
				register uint64_t U = tw[xx];
				U |= (uint64_t)twi[xx] << 56;
				t54g2.Add(U);
			}
		}
	}
	uint64_t* tdump = t54g2.t128[0].t;
	for (int i = 0; i < 10; i++)
		cout << Char54out(tdump[i]) << " " << (tdump[i] >> 56) << endl;
}
void T54G3::BuildG3() {//extract, shrink, sort g2
	Init();
	BF128 vsize[25][30];
	memset(vsize, 0, sizeof vsize);
	uint64_t tw[128 * 30];// can accept 40 uas per i81
	uint32_t twi[128 * 30], twn[128 * 30]; // i81,size 
	uint32_t ntw = 0;
	{ // extract still valid
		register uint64_t F = spb_0_15[7].all_previous_cells,
			A = spb_0_15[7].active_cells;
		cout << Char54out(F) << "F " << endl;
		cout << Char54out(A) << "A " << endl;
		for (uint32_t i81 = 0; i81 < 81; i81++) {
			GUAH54::GUA& gt = guah54.tg3[i81];
			if (gt.killer & F)continue; // killed or not valid
			uint32_t n1 = gt.nua;
			//cout << "i81=" << i81 << " n1=" << n1 << endl;
			if (n1 == 1) {// killer is enough
				register uint64_t U = gt.tua[0] & A;
				twi[ntw] = i81;
				twn[ntw] = (uint32_t)_popcnt64(U);
				tw[ntw++] = U;
				continue;
			}

			// must go in details and clean redundant
			if (n1 > 40) n1 = 40;
			uint64_t tw2[40];// valid
			uint32_t ntw2 = 0, iw2;
			uint64_t vsize2[25];// sorting tw2 by size
			memset(vsize2, 0, sizeof vsize2);

			// store still valid for active cells
			for (uint32_t i = 0; i < n1; i++) {
				register uint64_t U = gt.tua[i];
				if (U & F)continue; // hit
				U &= A;// keep active cells
				uint64_t cc = _popcnt64(U);
				vsize2[cc] |= (uint64_t)1 << ntw2;
				tw2[ntw2++] = U;
			}
			if (vsize2[0]) {// force true (one not hit)
				twi[ntw] = i81;
				twn[ntw] = 0;
				tw[ntw++] = 0;
				continue;
			}
			// take back per size and send to main table
			uint32_t istart = ntw; // start for redundancy
			for (int i = 0; i < 20; i++) {
				register uint64_t V = vsize2[i];
				while (bitscanforward64(iw2, V)) {
					V ^= (uint64_t)1 << iw2;
					register uint64_t U = tw2[iw2], nU = ~U;
					for (uint32_t j = istart; j < ntw; j++)
						if (!(tw[j] & nU)) { U = 0; break; }
					if (U) {
						twi[ntw] = i81;
						twn[ntw] = i;
						tw[ntw++] = U;
					}
				}
			}
		}// end i81 
	}// end build tw
	//cout << " buildg2 ntw=" << ntw << endl;


	// sort tw by size 
	for (uint32_t i = 0; i < ntw; i++) {
		uint32_t bloc = i >> 7, ir = i - (bloc << 7);
		vsize[twn[i]][bloc].setBit(ir);
	}
	// send the sorted table tw to t54g2
	cout << "BuildG3 back g3 ntw=" << ntw << endl;

	uint32_t nbl64 = (ntw + 63) >> 6, x, n = 0;

	for (int i1 = 0; i1 < 25; i1++) {
		uint64_t* tb64 = vsize[i1]->bf.u64;
		for (uint32_t i2 = 0; i2 < nbl64; i2++) if (tb64[i2]) {
			register uint64_t V = tb64[i2];
			while (bitscanforward64(x, V)) {
				V ^= (uint64_t)1 << x;
				uint32_t xx = x + (i2 << 6);
				register uint64_t U = tw[xx];
				U |= (uint64_t)twi[xx] << 56;
				t54g3.Add(U);
			}
		}
	}
	uint64_t* tdump = t54g3.t128[0].t;
	for (int i = 0; i < 10; i++)
		cout << Char54out(tdump[i]) << " " << (tdump[i] >> 56) << endl;
}


inline void G17B::BuildGua(BF128& w, int cc) {
	register uint64_t ua12 = w.bf.u64[0];
	uint64_t ua54 = (ua12 & BIT_SET_27) |
		((ua12 & BIT_SET_B2) >> 5);
	w.bf.u64[0] = ua54;// store in 54 cells mode
	if (cc < 4) {
		register uint32_t A = w.bf.u32[2];
		if (cc == 2) 	w.bf.u32[2] = genb12.bands3[0].GetI81_2(A);
		else w.bf.u32[2] = genb12.bands3[0].GetI81_3(A);
	}
}

inline int G17B::BuildGua(BF128& w) {// sockets2 sockets 3
	STD_B3& b = genb12.bands3[cur_ib];
	register uint32_t A = w.bf.u32[2];
	int cc=  _popcnt32(A);
	register uint64_t ua12 = w.bf.u64[0];
	uint64_t ua54 = (ua12 & BIT_SET_27) |
		((ua12 & BIT_SET_B2) >> 5);
	w.bf.u64[0] = ua54;// store in 54 cells mode
	if (cc > 6) return cc;// not one socket
	if (cc < 2) {
		w.bf.u32[2] = b.GetI81_2(A);return 2;}
	if (cc == 3) {
		w.bf.u32[2] = b.GetI81_3(A);return 2;	}
	// is it a socket 2 one stack 2 columns
	register uint32_t B = (A | (A >> 9) | (A >> 18)) & 0777;// columns
	int cc2 = _popcnt32(B);
	if (cc2 > 4) return cc;// not one socket
	if (cc2 == 4) {// can be one stack empy not one socket
		if( (!(B & 07007007)) || (!(B & 070070070)) ||
		   (!(B & 0700700700)) )return cc;	}
	int stack,cell1,cell2;// this is a gua 4/6
	if (_Popcount(B & 07007007) == 2) {	A &= 07007007; stack = 0;}
	else if (_Popcount(B & 070070070) == 2) { A &= 070070070; stack = 1; }
	else { A &= 0700700700; stack = 2; }
	bitscanforward(cell1, A);	bitscanreverse(cell2, A);
	int dig1 = b.band0[cell1], dig2 = b.band0[cell2],
		digs=(1<<dig1) | (1<<dig2),
		istart = 27 * stack, iend = istart + 27;
	for (int i = istart; i < iend; i++) {
		if (tsgua2[i].digs == digs) {
			w.bf.u32[2] = i;
			return 2;
		}
	}
	return cc;// should never be
}




/*
inline void GCHK::BuildGua(BF128 & w, int cc) {
	register uint64_t ua12 = w.bf.u64[0];
	uint64_t ua54 = (ua12 & BIT_SET_27) |
		((ua12 & BIT_SET_B2) >> 5);
	w.bf.u64[0] = ua54;// store in 54 cells mode
	if (cc < 4) {
		register uint32_t A = w.bf.u32[2], B, i;
		if (cc == 2) 	w.bf.u32[2] = GetI27(A);
		else {
			for (i = 0, B = 7; i < 9; i++, B <<= 3)
				if (A&B) { w.bf.u32[2] = i; break; }
		}
	}
}


//________ start expand bands 1+2 3 steps


void G17B::Init7p_guas() {
	register uint64_t Ac = spb_0_15[7].active_cells;
	chunkh.ApplyClean(tclues, 6);
	// _______________extract guas2 and vectors
	BF128 t[2000], tw;// temp storage for still active
	uint32_t  tt[15][100], ntt[15],
		nt = 0, n2 = 0, n3 = 0, nm = 0;
	memset(ntt, 0, sizeof ntt);
	// split guas2   in 128 vectors ( room for more )
	for (int i = 0; i < NG23_128; i++) gx_128x[i].Init();
	for (uint32_t ic2 = 0; ic2 <= chunkh.ic2; ic2++) {
		CHUNK3B& w = chunkh.c2[ic2];
		register uint64_t V = w.vclean;
		register uint32_t ir;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			register uint32_t i27 = w.tu27[ir];
			register uint64_t pat12 = w.tu54[ir] & Ac;
			tw.bf.u64[0] = pat12;
			tw.bf.u32[3] = i27;
			tw.bf.u32[2] = 1 << i27;
			register uint32_t cc = (uint32_t)_popcnt64(pat12);
			// below cc=10 load direct else sort by size
			if (cc < 10) { g_128h.g2_128[0].Add(tw); n2++; continue; }
			if (cc < 20) cc -= 10;
			else if (nt < 500)cc = 9;
			else continue;
			tt[cc][ntt[cc]++] = nt;
			if (nt < 2000)t[nt++] = tw;
			else break;
		}
	}


	g_128h.GetStart(n2, n3, nm);

}

/*


	for (uint32_t i = 0; i < 10; i++) {
		uint32_t* tti = tt[i];
		for (uint32_t j = 0; j < ntt[i]; j++) {
			BF128 w = t[tti[j]];
			if (n2++ < 256)g2_256[0].Add(w);
			else g2_256[1].Add(w);
		}
	}
	// ___________________extract guas3
	for (uint32_t ic3 = 0; ic3 <= chunkh.ic3; ic3++) {
		CHUNK3B& w = chunkh.c3[ic3];
		register uint64_t V = w.vclean;
		register uint32_t ir;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			uint32_t i9 = w.tu27[ir];
			register uint64_t pat12 = w.tu54[ir] & Ac;
			tw.bf.u64[0] = pat12;
			tw.bf.u32[3] = i9;
			tw.bf.u32[2] = 1 << i9;
			g3_256.Add(tw);
		}
	}
	tw.bf.u32[3] = 0;
	// _____________________________________________extract guas4
	for (uint32_t ic4 = 0; ic4 <= chunkh.ic4; ic4++) {
		CHUNK3B& w = chunkh.c4[ic4];
		register uint64_t V = w.vclean;
		register uint32_t ir;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			tw.bf.u64[0] = w.tu54[ir] & Ac;
			tw.bf.u32[2] = w.tu27[ir];
			gm_256[0].Add(tw);
			nm++;
		}
	}
	// _____________________________________________extract guas5
	for (uint32_t ic5 = 0; ic5 <= chunkh.ic5; ic5++) {
		CHUNK3B& w = chunkh.c5[ic5];
		register uint64_t V = w.vclean;
		register uint32_t ir;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			tw.bf.u64[0] = w.tu54[ir] & Ac;
			tw.bf.u32[2] = w.tu27[ir];
			gm_256[0].Add(tw);
			nm++;
		}
	}
	// _____________________________________________extract guasmore
	for (uint32_t icm = 0; icm <= chunkh.icmore; icm++) {
		CHUNK3B& w = chunkh.cmore[icm];
		register uint64_t V = w.vclean;
		register uint32_t ir;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			tw.bf.u64[0] = w.tu54[ir] & Ac;
			tw.bf.u32[2] = w.tu27[ir];
			if (nm++ < 256)gm_256[0].Add(tw);
			else gm_256[1].Add(tw);
		}
	}
	g_128h.GetStart(n2, nm);
}
*/

/*
uint32_t tbuilexpand3[12][500];// 6+6

int GCHK::BuildB2Table() {// extract and reorder still possible guas
	uint32_t is1 = 0, ntt[6], i27;
	memset(ntt, 0, sizeof ntt);
	register uint32_t AC = scritb3.active,
		F = scritb3.assigned;
	{
		register uint32_t V = scritb3.pairs27;
		while (bitscanforward64(i27, V)) {
			V ^= 1 << i27;
			tbuilexpand3[2][ntt[2]++] = tg2[i27].pat;
		}
		V = scritb3.minix[0];// triplets
		if (V) {// add guas3 if any (usually 0)
			for (int i = 0; i < 9; i++)if (V & (1 << i))
				tbuilexpand3[3][ntt[3]++] = 7 << (3 * i);
		}

	}

	register uint32_t * t2 = tbuilexpand3[2];
	for (uint32_t i = 0; i < ntw3; i++) {
		register uint32_t U = tw3[i];
		if (!(U&F)) {
			U &= AC;
			int ccu = _popcnt32(U);
			if (!ccu) return 1;
			if (ccu == 1) { is1 |= U; continue; }
			{// clear redundancy
				register uint32_t nu = ~U, x = 0;
				for (uint32_t j = 0; j < ntt[2]; j++)
					if (!(nu&t2[j])) { x = 1; break; }
				if (x) continue;
			}
			if (ccu > 5) ccu = 5;
			tbuilexpand3[ccu][ntt[ccu]++] = U;
		}
	}
	// add small band 3 residual uas
	for (uint32_t i = 0; i < myband3.nua; i++) {
		register uint32_t U = myband3.tua[i];
		if (!(U&F)) {
			U &= AC;
			int ccu = _popcnt32(U);
			if (!ccu) return 1;
			if (ccu == 1) { is1 |= U; continue; }
			{// clear redundancy
				register uint32_t nu = ~U, x = 0;
				for (uint32_t j = 0; j < ntt[2]; j++)
					if (!(nu&t2[j])) { x = 1; break; }
				if (x) continue;
			}
			if (ccu > 5) ccu = 5;
			tbuilexpand3[ccu][ntt[ccu]++] = U;
		}
	}


	for (uint32_t i = 0; i < ntaddgob3; i++) {
		register uint32_t U = taddgob3[i];
		if (!(U&F)) {
			U &= AC;
			int ccu = _popcnt32(U);
			if (!ccu) return 1;
			if (ccu == 1) { is1 |= U; continue; }
			{// clear redundancy
				register uint32_t nu = ~U, x = 0;
				for (uint32_t j = 0; j < ntt[2]; j++)
					if (!(nu&t2[j])) { x = 1; break; }
				if (x) continue;
			}
			if (ccu > 5) ccu = 5;
			tbuilexpand3[ccu][ntt[ccu]++] = U;
		}
	}
	ntw3_2 = is1_tw3=0;
	if (is1) {// now singles to assign
		if (scritb3.AddAssign(is1)) 	return 2;
		AC = scritb3.active;	F = scritb3.assigned;
		int cc = _popcnt32(F);		if (cc > ncluesb3)return 1;
		// after added clues, revise the others
		uint32_t  ntt2[6];//tt2[i][500]is wubuf.tbuilexpand3[i+5][500]
		memset(ntt2, 0, sizeof ntt2);
		for (int i = 2; i < 6; i++)
			for (uint32_t j = 0; j < ntt[i]; j++) {
				register uint32_t U = tbuilexpand3[i][j];
				if (!(U&F)) {
					U &= AC;
					register int cc = _popcnt32(U);
					if (!cc) return 3; // dead branch
					if (cc > 5) cc = 5;
					tbuilexpand3[cc + 6][ntt2[cc]++] = U;
				}
			}
		if (ntt2[1]) {// send back the bit field status
			for (uint32_t j = 0; j < ntt2[1]; j++)
				is1_tw3 |= tbuilexpand3[7][j];
		}
		for (int i = 1; i < 6; i++)
			for (uint32_t j = 0; j < ntt2[i]; j++) {
				// clear redundancy
				register uint32_t U = tbuilexpand3[i + 6][j],
					nu = ~U, x = 0;
				for (uint32_t i = 0; i < ntw3_2; i++)
					if (!(nu&tw3_2[i])) { x = 1; break; }
				if (x) continue;
				tw3_2[ntw3_2++] = U;
			}
	}
	else {// apply directly stored uas
		for (int i = 2; i < 6; i++)
			for (uint32_t j = 0; j < ntt[i]; j++) {
				// clear redundancy
				register uint32_t U = tbuilexpand3[i][j],
					nu = ~U, x = 0;
				for (uint32_t k = 0; k < ntw3_2; k++)
					if (!(nu & tw3_2[k])) { x = 1; break; }
				if (x) continue;
				tw3_2[ntw3_2++] = U;
			}
	}

	return 0;
}
*/