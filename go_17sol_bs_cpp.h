#ifdef XXXXXX
//___ start process expand bands collect uas guas ...
const char * diagband =  "364912578715834962892675341"; // la bonne
const char * diagband3 = "276341895538297416941568237";
const char * diagpuz = "...45..........6...89...1....4......7.5.........6..3.......1.95.3.2.............7";
#endif

#define TEST_ON 1
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
					if (Bf & 0777) { nr++; bfr |= Bf & 0777; irow = 0; }// row1
					if (Bf & 0777000) { nr++; bfr |= Bf & 0777000; irow = 1; }
					if (Bf & 0777000000) { nr++; bfr |= Bf & 0777000000; irow = 2; }
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
						b3.g.ua2bit27[i] = bit27;
						b3.i_27_to_81[i27] = i;
						b3.i_81_to_27[i] = i27;
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
	cout << " grid to use in ua collector" << endl;

	//StartInitDebug();
	//for (int i = 0; i < genb12.nband3; i++) {
		//STD_B3& b3 = genb12.bands3[i];
		//cout << b3.band << " b3 i=" << i << endl;
		//for (int i = 0; i < (int)b3.nua; i++)
		//	cout << Char27out(b3.tua[i]) << endl;
	//}

}

void G17B::Start() {// processing an entry 
	if (aigstop)return;
#ifdef TEST_ONEB2
	if (p_cpt2g[0]) return;
#endif
	is_test_on = 0;
	p_cpt2g[0] ++;
	p_cpt2g[1] += genb12.nband3;
	if (p_cpt2g[0] > 1) {		aigstop = 1; return;	}
	// back if min 66 or not 5 clues in band 2  pass 3

	StartInit();// do all preliminary setups
#ifdef TEST_ON
	StartPrint();
#endif
	UaCollector();//_do uas/guas2_3 initial harvest
	//return;
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
struct EXTUAS {// get a sorted subset of uas for next "collect" step
	uint64_t  t2[300];
	uint32_t  nt2;
	void GetSortB(uint64_t* t, uint32_t n, uint64_t filter, uint32_t nold = 0) {
		BF128 vsize[25][4];
		uint64_t  t2a[512];
		uint32_t  nt2a = 0;
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

struct GUAH {// handler for initial collection of guas 2 3

	struct GUA {
		uint64_t tua[128];
		uint32_t nua, type, i81;
		inline void Add(uint64_t ua) {
			if (nua < 128)tua[nua++] = ua;
		}
		inline void AddIf(uint64_t ua) {
			for (uint32_t i = 0; i < nua; i++) {
				if (ua == tua[i]) return;
			}
			tua[nua++] = ua;
		}

		inline void Init(uint32_t n, uint32_t t, uint32_t i) {
			nua = n; type = t; i81 = i;
		}
		int Load(uint64_t* tu, uint64_t bf) {
			int n = 0;
			register uint64_t nbf = ~bf;
			for (uint32_t i = 0; i < nua; i++) {
				register uint64_t U = tua[i], cc = _popcnt64(U);
				if (n > 10) break;// 
				if (n > 5 && cc > 12) break;// 
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
				register uint64_t  cc = _popcnt64(tua[i] & BIT_SET_2X);
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
				cout << Char2Xout(tua[i]) << " " << _popcnt64(tua[i] & BIT_SET_2X)
					<< " " << i << endl;
			}
		}
	}tg2[81], tg3[81], guaw;
	void Init() {
		for (int i = 0; i < 81; i++) {
			tg2[i].Init(0, 0, i); tg3[i].Init(0, 1, i);
		}
	}
	inline void Add2(uint64_t u, uint32_t i) { if (_popcnt64(u) < 18)		tg2[i].Add(u); }
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
		for (int i = 0; i < 81; i++) if (tg2[i].nua) {
			tg2[i].Debug(0);
		}
	}
	void Dump2all3() {
		for (int i = 0; i < 81; i++) if (tg3[i].nua) {
			tg3[i].Debug(0);
		}
	}
}guah;

void G17B::UaCollector() {
	zh2gxn.SetupFsol(grid0);
	zh2b[0].InitBands12(grid0);
	zh2gxn.InitKnown(extuas.t2, &extuas.nt2);
	FirstUasCollect();// uas 2 digits bands and stacks
	//return;
//	aigstop = 1; return;
	SecondUasCollect();// uas 345 bands 1+2
	UasCollect6_7();
	Guas2Collect();
	UasCollect4box();// 4 box in bands 1+2
	zh2b[0].InitBands12(grid0);// restore 
	//________________________________ insert stack known uas 
	STD_B3 bw;
	BANDMINLEX::PERM pp;
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		STD_B3& mb3 = genb12.bands3[ib3];
		memcpy(&grid0[54], mb3.band0, sizeof genb12.bands3[0].band0);
		int grid_diag[81];
		char ze[81];
		for (int i = 0; i < 81; i++)grid_diag[C_transpose_d[i]] = grid0[i];
		for (int i = 0; i < 81; i++)  ze[i] = grid_diag[i] + '1';
		for (int istack = 0; istack < 3; istack++) {
			int ir = bandminlex.Getmin(&grid_diag[27 * istack], &pp, 0);
			if (ir < 0) {//would be bug} 
				cout << "look for stack bad return" << endl;
				return;
			}
			bw.InitBand2_3(pp.i416, &ze[27 * istack], pp, istack);
			//bw.PrintStatus();
			for (uint32_t i = 0; i < bw.nua; i++) {
				BF128 wd; wd.SetAll_0();// 3x mode to build guam
				uint32_t  zed = 0, ua = bw.tua[i]&BIT_SET_27,ccua=_popcnt32(ua);
				for (uint32_t j = 0; j < 27; j++) {
					if (ua & (1 << j)) {
						int cell = C_transpose_d[j + 27 * istack];
						wd.Set_c(cell);
					}
				}
				int cc3 = _Popcount(wd.bf.u32[2]);
				//cout << Char27out(zed) << endl;
				//cout << Char2Xout(wd.bf.u64[0]) << " ";
				//cout << Char27out(wd.bf.u32[2])<<" "<<wd.Count()<<" "<<cc3 << endl;
				if (cc3 > 3)mb3.Addguam(wd);
				else if(cc3==2  &&  ccua>10) { // must build g2
					int i81=mb3.GetI81_2(wd.bf.u32[2]);
					guah.tg2[i81].AddIf(wd.bf.u64[0]);
				}
				else if (cc3 == 3  &&  ccua>9) { // must build g3
					int i81 = mb3.GetI81_3(wd.bf.u32[2]);
					guah.tg3[i81].AddIf(wd.bf.u64[0]);
				}

			}
		}
		//cout << "band summary ib3 =" << ib3 << "\t";
		//mb3.DumpGuam();

	}
	//zh2b_g.Init_g0(grid0);// reinit brute force bands 1+2
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
				if (cct <= cc)continue; 
				// insert here and check super set
				for (uint32_t jua = ntall; jua > iua; jua--)
					tall[jua] = tall[jua - 1];
				tall[iua] = w;// new inserted
				tall[iua].bf.u32[3] = floor | (cc << 16);
				ntall++;
				return 2;
			}
			w.bf.u32[3] = floor | (cc << 16);
			tall[ntall++] = w;// new added
			return 1;
		}
	}tuas81;

	cout << "entry FirstUasCollect()" << endl;
	for (int i = 0; i < 36; i++) {
		uint32_t myf = floors_2d[i];
		zh2_2[0].GoZ2D(myf);
		//continue;
		if (zh2gxn.nua) {
			//cout << "found UAs for i=" << i << endl;
			for (uint32_t i = 0; i < zh2gxn.nua; i++) {
				register uint64_t w = zh2gxn.tua[i], cc = _popcnt64(w);
				//cout << Char2Xout(w) << " " << cc << endl;
			if (cc > 15)cc = 15;
			w |= cc << 59;
			tuasb12.AddInit(w, myf, 2);
			}
		}

	}

	// do zhou2[0].GoZ2(myf) for al bands looing for guas
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
		STD_B3& mb3 = genb12.bands3[ib3];
		memcpy(&grid0[54], mb3.band0, sizeof genb12.bands3[0].band0);
		zhgxn.SetupFsol(grid0);
		//cout << "ib3 =" << ib3 << endl;
		for (int i = 0; i < 36; i++) {
			uint32_t myf = floors_2d[i];
			zhou2[0].GoZ2(myf);
			//cout << "back nua=" << zhgxn.nua << endl;
			if (zhgxn.nua == 1) {// this is a 18 cells 12+6
				BF128 w = zhgxn.tua[0];
				register uint32_t C = w.bf.u32[2];
				uint32_t  ncol = _popcnt32((C | (C >> 9) | (C >> 18)) & 0777);
				if (ncol != 4)mb3.Addguam(w);
				else {// this is a 6 clues gua2
					//cout << Char2Xout(w.bf.u64[0]) << "socket 2 12+6" << endl;
					//for (int i = 0; i < 3; i++) {
					//	cout << Char9out(C) << endl;
					//	C >>= 9;
					//}
				}
			}
			else {
				for (uint32_t j = 0; j < zhgxn.nua; j++) {
					BF128 w = zhgxn.tua[j];
					int cc = w.Count();
					int aig = 1;// check for subsets
					if (cc == 18)continue; // has a subset
					for (uint32_t k = 0; k < zhgxn.nua; k++) if (k != j) {
						BF128 w2 = zhgxn.tua[k];
						if (w2.isSubsetOf(w)) { aig = 0; break; }
					}
					if (aig) {
						uint64_t w1 = w.bf.u64[0], cc = _popcnt64(w1);
						register uint32_t C = w.bf.u32[2];
						uint32_t cc2 = _popcnt32(C),
							ncol = _popcnt32((C | (C >> 9) | (C >> 18)) & 0777),
							nplug = 2 * ncol - cc2;
						if ((!cc) || (!cc2)) continue;//ua bands 1+2 or ua band 3
						if (nplug > 2) { mb3.Addguam(w); continue; }
						//cout << Char2Xout(w.bf.u64[0]) << "socket 2 "<< cc<<" + "<<cc2 << endl;
						//for (int i = 0; i < 3; i++) {
						//	cout << Char9out(C) << endl;
						//	C >>= 9;
						//}
					}
				}
			}		
		}
		//mb3.DumpGuam2d();
	}



}
void G17B::SecondUasCollect() {// collect 345 digits in bands 1+2
	//cout << "seconuacollect" << endl;
	//zh2gxn.InitKnown(tuasb12.t2, &tuasb12.nt2);
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
		//if (extuas.nt2 && HasUaHitFalse(extuas.t2,0, extuas.nt2, i81))
		//	continue;

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
			//if (extuas.nt2>istart && HasUaHitFalse(extuas.t2, istart,extuas.nt2, i81))
			//	continue;
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
			//if (extuas.nt2 > istart && HasUaHitFalse(extuas.t2, istart, extuas.nt2, i81))
			//	continue;
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
			//if (diag2 == 2)extuas.Dump();
			zh2gxn.nkguas = extuas.nt2;
			uint32_t istart = extuas.nt2;
			extuas.nt2 += ttt1.Load(&extuas.t2[istart], 3, i);
			//if (diag2 == 2)extuas.Dump();
			// check no ua hitting cells "false"
			//if (extuas.nt2 > istart && HasUaHitFalse(extuas.t2, istart, extuas.nt2, i81))
			//	continue;
			int ir = zh2_5[0].GoZ5G2(fl, w.col1, w.dig1, w.col2, w.dig2), ir2;
			if (ir < 0)continue; // locked
			//if (diag2 == 2)cout << "back goz5g2 ir=" << ir << endl;
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
/*
    seems very expensive 
	//___ find guas2 6 digits
	for (int i81 = 0; i81 < 81; i81++) if (gsock2.On(i81)) {
		//int diag = 0; if (i81 == 6)diag = 1;
		if (guah.IsUa4(i81)) continue;// nothing to do
		SGUA2 w = tsgua2[i81];
		if (diag == 1)cout << Char9out(w.digs) << "digs missed gua2 bug analysis" << endl;
		GUAH::GUA& gt = guah.tg2[i81];
		if (gt.nua >= 30) gt.nua = 30;// keep room for short more uas
		for (int i = 0; i < 84; i++) {
			int fl = 0777 ^ floors_3d[i];
			if (!((fl & w.digs) == w.digs)) continue;
			//int diag2 = 0;
			//if (diag == 1 && fl == 0363) diag2 = 2;//125678
			//if (diag2 == 2) cout << Char9out(fl) << "missed gua2 bug analysis" << endl;

			uint64_t isfl = 0;// cells ok of the floors
			for (int i2 = 0, bit = 1; i2 < 9; i2++, bit <<= 1)
				if (fl & bit)	isfl |= zh2gxn.fsol[i2];
			extuas.nt2 = gt.Load(extuas.t2, isfl);
			zh2gxn.nkguas = extuas.nt2;
			uint32_t istart = extuas.nt2;
			extuas.nt2 += ttt1.Load(&extuas.t2[istart], 3, i);
			//if (diag2 == 2)extuas.Dump();
			// check no ua hitting cells "false"
			//if (extuas.nt2 > istart && HasUaHitFalse(extuas.t2, istart, extuas.nt2, i81))
			//	continue;
			zh2b[0].InitBfG2(BIT_SET_2X^ isfl, w.col1, w.dig1, w.col2, w.dig2);
			//continue;
			zh2b[0].Dob(14);
			//if (diag2 == 2)cout << "back goz5g2 zh2gxn.nua=" << zh2gxn.nua << endl;
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

	*/
	ng2=guah.CutG2(30);
	if (is_test_on)cout << "ng2 after cut30=" << ng2 << endl;
	//guah.Dump2all2();

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

	// here insert all missing uas one band (1 or 2)
	for (int i = 0; i < 2; i++) {
		STD_B416& wb = (i) ? myband2 : myband1;
		BF64 wt; wt.bf.u64 = 0;
		uint32_t* tu = wb.tua, nu = wb.nua;
		for (uint32_t j = 0; j < nu; j++) {
			register uint32_t u = tu[j] & BIT_SET_27, cc = _popcnt32(u);
			wt.bf.u32[i] = u;
			if (cc > 12)if(n< UA12SIZE)t[n++] = wt.bf.u64;
		}
	}	
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
							//cout << Char54out(R) << " erased" << endl;;
							//cout << Char54out(tw[i]) <<" i="<<i << endl;
							R = 0; break;						}
					}
					if (R)tw[ntw++] = R;
					else continue;
				}
				AddA(R);
			}
		}
	}
	if (g17b.is_test_on > 1)cout << "ntw=" << ntw << " nold" << n << endl;
	//if(ntw!=n)ta128[0].Dump();
}
int T54B12::Build_tb128() {
	BF128 tvw[20];
	uint32_t lastbloc = t54b12.nablocs;
	tvw[0] = spb_0_15[3].v;
	for (uint32_t i = 1; i <= lastbloc; i++) {
		TUVECT& vv = ta128[i];
		BF128 v = vv.v0, * vc = vv.vc;
		for (uint32_t ic = 0; ic < 3; ic++)
			v &= vc[g17b.tclues[ic]];
		tvw[i] = v;
	}

	// apply active on still valid uas and flag by size
	BF128 vsize[19][UABNBLOCS];
	memset(vsize, 0, sizeof vsize);
	//uint32_t nbl64 = (n + 63) >> 6, x;
	uint64_t tw[128 * UABNBLOCS];
	uint32_t ntw = 0;
	{
		register uint64_t Ac = spb_0_15[3].active_cells;
		for (uint32_t i = 0; i <= lastbloc; i++) {
			register uint64_t* t = ta128[i].t;
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
	if (g17b.is_test_on > 1)cout << "ua3c ntw=" << ntw << " nb128=" << nb128 << endl;

	return 0;
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
				if (t54b12.IsNotRedundant(U))
					AddC(U);
			}
		}
	}
	if (g17b.is_test_on > 1)cout << "ua6c ntw=" << ntw << " nc128=" << nc128 << endl;
	return 0;
}


void GUAH54::Build() {// cut to 30 switch to 54 killer

	uint64_t* pbuf = gbuf;
	for (int i81 = 0; i81 < 81; i81++) {
		tg2[i81].Init(pbuf, 0, i81);
		if (g17b.gsock2.On(i81)) {
			GUAH::GUA& g0 = guah.tg2[i81];
			GUA54& gd = tg2[i81];
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
			GUA54& gd = tg3[i81];
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
void GUAH54::Build2(uint64_t filter, uint64_t active) {
	register uint64_t F = filter,
		A = active;

	uint64_t* pbuf = gbuf;
	for (int i81 = 0; i81 < 81; i81++) {
		tg2[i81].Init(pbuf, 0, i81);
		if (g17b.gsock2.On(i81)) {
			GUA54& g0 = guah54.tg2[i81];
			GUA54& gd = tg2[i81];
			if (g0.killer & F)continue; // killed or not valid
			uint32_t n1 = g0.nua;
			if (n1 == 1) {// killer is enough
				register uint64_t U = g0.tua[0] & A;
				gd.Add(U);
				pbuf++;
				continue;
			}
			uint64_t vsize[25];// sorting tw2 by size
			memset(vsize, 0, sizeof vsize);
			uint64_t tw[40];// valid
			uint32_t ntw = 0, iw;
			// store still valid for active cells
			for (uint32_t i = 0; i < n1; i++) {
				register uint64_t U = g0.tua[i];
				if (U & F)continue; // hit
				U &= A;// keep active cells
				uint64_t cc = _popcnt64(U);
				vsize[cc] |= (uint64_t)1 << ntw;
				tw[ntw++] = U;
				if (ntw >= 40) break;
			}
			if (vsize[0]) {// force true (one not hit)
				//cout << "0 g2 found i81=" << i81 << endl;
				gd.Add(0);
				pbuf++;
				continue;
			}			// take back per size and send to main table
			uint32_t istart = ntw; // start for redundancy
			for (int i = 0; i < 20; i++) {
				register uint64_t V = vsize[i];
				while (bitscanforward64(iw, V)) {
					V ^= (uint64_t)1 << iw;
					gd.AddCheck(tw[iw]);
				}
			}
			gd.nuamax = gd.nua + 5;
			pbuf += gd.nuamax;// lock the storing place
		}
	}
	for (int i81 = 0; i81 < 81; i81++) {
		tg3[i81].Init(pbuf, 1, i81);
		if (g17b.gsock3.On(i81)) {
			GUA54& g0 = guah54.tg3[i81];
			GUA54& gd = tg3[i81];
			if (g0.killer & F)continue; // killed or not valid
			uint32_t n1 = g0.nua;
			if (n1 == 1) {// killer is enough
				register uint64_t U = g0.tua[0] & A;
				gd.Add(U);
				pbuf++;
				continue;
			}
			uint64_t vsize[25];// sorting tw2 by size
			memset(vsize, 0, sizeof vsize);
			uint64_t tw[30];// valid
			uint32_t ntw = 0, iw;
			// store still valid for active cells
			for (uint32_t i = 0; i < n1; i++) {
				register uint64_t U = g0.tua[i];
				if (U & F)continue; // hit
				U &= A;// keep active cells
				uint64_t cc = _popcnt64(U);
				vsize[cc] |= (uint64_t)1 << ntw;
				tw[ntw++] = U;
				if (ntw >= 30) break;
			}
			if (vsize[0]) {// force true (one not hit)
				//cout << "0 g3 found i81=" << i81 << endl;
				gd.Add(0);
				pbuf++;
				continue;
			}			// take back per size and send to main table
			uint32_t istart = ntw; // start for redundancy
			for (int i = 0; i < 20; i++) {
				register uint64_t V = vsize[i];
				while (bitscanforward64(iw, V)) {
					V ^= (uint64_t)1 << iw;
					gd.AddCheck(tw[iw]);
				}
			}
			gd.nuamax = gd.nua + 5;
			pbuf += gd.nuamax;// lock the storing place
		}
	}

	//Dumpall2();
	//Dumpall3();
}
BF128 GUAH54::GetG2(uint64_t bf) {
	register uint64_t F = bf;
	BF128 g2;
	g2.SetAll_0();
	for (int i81 = 0; i81 < 81; i81++) {
		GUA54& g0 = tg2[i81];
		if (g0.killer & F)continue; // killed or not valid
		uint32_t n1 = g0.nua;
		if (n1 == 1) {	g2.setBit(i81);	continue;}
		for (uint32_t i = 0; i < n1; i++) {
			register uint64_t U = g0.tua[i];
			if (U & F)continue; // hit
			g2.setBit(i81);	// first not hit is ok
			break;
		}
	}
	return g2;
}
BF128 GUAH54::GetG3(uint64_t bf) {
	register uint64_t F = bf;
	BF128 g3;
	g3.SetAll_0();
	for (int i81 = 0; i81 < 81; i81++) {
		GUA54& g0 = tg3[i81];
		if (g0.killer & F)continue; // killed or not valid
		uint32_t n1 = g0.nua;
		if (n1 == 1) { g3.setBit(i81);	continue; }
		for (uint32_t i = 0; i < n1; i++) {
			register uint64_t U = g0.tua[i];
			if (U & F)continue; // hit
			g3.setBit(i81);	// first not hit is ok
			break;
		}
	}
	return g3;
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
		//if (_popcnt64(ac) < 45) {
			//cout << Char54out(sn->all_previous_cells) << " 3clues [3] " << p_cpt2g[3] << endl;
			if (t54b12.Build_tb128()) goto next;
			Expand_46();
		//}
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
		//if ( _popcnt64(ac) < 27) {
		//if (p_cpt2g[4]<3 ||_popcnt64(ac) < 27) {
			//cout << Char54out(sn->all_previous_cells) << "\t 6clues [4] " << p_cpt2g[4] 
			//	<<" "<< _popcnt64(ac) << endl;
			if (t54b12.Build_tc128()) goto next;
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
			if (mincluesb3 == 6)Expand_7_11();
			else 		GoExpand_7_10();
#endif
		//}
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
			GUA54  & gt = guah54.tg2[i81];
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
	// and setup vcl for next steps 
	for (int i = 0; i < NG2BLOCS6; i++) vcl[0][i] = t128[i].v0;
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
			GUA54 & gt = guah54.tg3[i81];
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
	// and setup vcl for next steps 
	for (int i = 0; i < NG3BLOCS6; i++) vcl[0][i] = t128[i].v0;
	uint64_t* tdump = t54g3.t128[0].t;
	for (int i = 0; i < 10; i++)
		cout << Char54out(tdump[i]) << " " << (tdump[i] >> 56) << endl;
}


/*
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
*/


int G17B::IsValid7p(SPB03* sn) {
	if (zh2b[1].IsValid(tclues, sn->ncl)) {
		anduab12 = ~0;
		register uint64_t ua54;
		for (uint32_t i = 0; i < zh2gxn.nua; i++) {
			register uint64_t ua = zh2gxn.tua[i], cc = _popcnt64(ua);
			ua54 = (ua & BIT_SET_27) | ((ua & BIT_SET_B2) >> 5);
			anduab12 &= ua54;
			if (cc < 21) {
				//tuasb12.tv128h.AddA(ua54);
				//tuasb12.tv128h.AddB(ua54);
				//tuasb12.tv128h.AddC(ua54);
			}
		}
		ua_ret7p = ua54;// return last (smaller)
		return 1;
	}
	return 0;
}
inline int G17B::GetNextCell(SPB03* s) {
	SPB03* sn;
	register int cell;
	uint64_t p = s->possible_cells;
	if (!p)return 1;
	bitscanforward64(cell, p);
	register uint64_t bit = (uint64_t)1 << cell;
	s->possible_cells ^= bit;
	tclues[s->ncl] = cell;
	s->active_cells ^= bit;
	sn = s + 1; *sn = *s; sn->ncl++;
	sn->cbs.Add(cell);
	sn->all_previous_cells |= bit;
	sn->v &= t54b12.tc128[0].vc[cell];
	return 0;
}
inline void G17B::GetNextUa(SPB03* sn) {
	register uint64_t  V;
	if ((V = sn->v.bf.u64[0])) {// next ua
		register uint32_t ir;
		bitscanforward64(ir, V);//relative index first active
		ua_ret7p = t54b12.tc128[0].t[ir];
	}
	else {// next ua must be here
		V = sn->v.bf.u64[1];
		register uint32_t ir;
		bitscanforward64(ir, V);//relative index first active
		ua_ret7p = t54b12.tc128[0].t[ir + 64];
	}
}
inline void G17B::GetNextUaAdd(SPB03* sn) {
	if (t54b12.nc128 <= 128) return;
	// more uas to check
	for (uint32_t i = 1; i <= t54b12.ncblocs; i++) {
		T54B12::TUVECT& vv = t54b12.tc128[i];
		BF128 v = vv.v0, * vc = vv.vc;
		for (uint32_t ic = 6; ic < sn->ncl; ic++)
			v &= vc[tclues[ic]];
		if (v.isNotEmpty()) {
			int ir2 = v.getFirst128();
			ua_ret7p = vv.t[ir2];
			return;
		}
	}
}
inline int G17B::GetLastAndUa(SPB03* sn, int d) {
	int aig = 0;
	register uint64_t  V, And = ~0;
	register uint32_t ir;
	if ((V = sn->v.bf.u64[0])) {
		aig = 1;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			And &= t54b12.tc128[0].t[ir];
		}
	}
	if (!And) return 1;
	if ((V = sn->v.bf.u64[1])) {
		aig = 1;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			And &= t54b12.tc128[0].t[ir + 64];
			if (!And) return 1;
		}
	}
	if (!And) return 1;
	if (t54b12.ncblocs) {
		// more uas to check
		for (uint32_t i = 1; i <= t54b12.ncblocs; i++) {
			T54B12::TUVECT& vv = t54b12.tc128[i];
			BF128 v = vv.v0, * vc = vv.vc;
			for (uint32_t ic = 6; ic < sn->ncl; ic++)
				v &= vc[tclues[ic]];
			if (v.isNotEmpty()) {
				aig = 1;
				int ir2;
				while ((ir2 = v.getFirst128()) >= 0) {
					v.clearBit(ir2);
					And &= vv.t[ir2];
				}
			}
		}
	}
	if (aig) ua_ret7p = And;
	return aig;
}



void DumpPotential(int det=0) {
	cout <<  "DumpPotential() " ;
	for (int i = 0; i < 6; i++) cout << ntbelow[i] << " ";
	cout << endl;
	if (!det) return;
	for (int i = 0; i < 6; i++)if (ntbelow[i]) {
		cout << Char54out(tandbelow[i]) << " and i" << i 
			<<" "<<_popcnt64(tandbelow[i]) << endl;
		if (i < 5) {
			if (det>1) {
				BF128* tb[5] = { tbelow7 ,tbelow8,tbelow9,tbelow10, tbelow11 };
				cout << "tbelow nclues=<i+7 ncl=" << i + 7 << endl;
				BF128* tbw = tb[i];
				for (uint32_t j = 0; j < ntbelow[i]; j++)
					cout << Char54out(tbw[j].bf.u64[0]) << " " << j << endl;
			}
			continue;
		}
		//for (uint32_t j = 0; j < ntbelow[i]; j++)
			//cout << Char54out(tfull[j]) << " " << j << endl;
		return;
	}

}

int G17B::Expand_7_10() {
	memset(ntbelow, 0, sizeof ntbelow);//7 8 9 10 11 full
	memset(tandbelow, 255, sizeof tandbelow);//7 8 9 10 11 full
	uint32_t nfull = 0, nbelow = 0;
	SPB03* sl = &spb_0_15[8], * s = sl, * sn;
	T54B12::TUVECT& tuv128 = t54b12.tc128[0];
	uint64_t* twu = tuv128.t;
	if (t54b12.nc128 < 128)t54b12.nc128 = 128;// force adds outside
	*s = spb_0_15[7];	// duplicate 6 for new vector
	s->possible_cells = twu[0];
	s->v = tuv128.v0;// initial nothing done
	//_______ start search clues 7_10
next:	// catch and apply cell in bitfields
	if (GetNextCell(s))	if (--s >= sl)goto next;	else {
		return (nfull | nbelow);
	}
	sn = s + 1;
#ifdef HAVEKNOWN
	if (!((~pk54) & sn->all_previous_cells)) {
		cout << Char54out(sn->all_previous_cells) << " on the path 7_11 ncl=" << sn->ncl << endl;
	}
#endif
	ua_ret7p = 0;
	if (sn->ncl == 10) {// 10 cells
		if (sn->cbs.IsFilt10_17()) goto next;
		GetNextUaAdd(sn);// check adds
		if (ua_ret7p)goto next;
		clean_valid_done = 0;
		p_cpt2g[5]++;
		tfullcbs[ntbelow[5]] = sn->cbs;
		tfull[ntbelow[5]++] = sn->all_previous_cells;
		tandbelow[5]&= sn->all_previous_cells;
		nfull++;
		goto next;
	}
	if (sn->ncl == 9) {// last step
		if (GetLastAndUa(sn)) // last not empty
			if (!ua_ret7p) goto next;
	}
	else {
		if (sn->v.isNotEmpty())GetNextUa(sn);// first 128
		else GetNextUaAdd(sn);
	}
	if (!ua_ret7p) {// add ua or go below
		if (IsValid7p(sn)) {// got uas to use
			if (sn->ncl == 10)ua_ret7p = anduab12;
		}
		else {// valid below
			clean_valid_done = 1;
			int isize = sn->ncl - 7;//table index 0 7 clues
			p_cpt2g[54+isize]++;
			BF128 w;
			w.bf.u64[0] = sn->all_previous_cells;
			w.bf.u64[1] = sn->active_cells;
			switch (isize) {
			case 0: tbelow7cbs[ntbelow[0]] = sn->cbs;  tbelow7[ntbelow[0]++] = w;
				tandbelow[0] &= sn->all_previous_cells; break;
			case 1: tbelow8cbs[ntbelow[1]] = sn->cbs; tbelow8[ntbelow[1]++] = w;
				tandbelow[1] &= sn->all_previous_cells; break;
			case 2: tbelow9cbs[ntbelow[2]] = sn->cbs; tbelow9[ntbelow[2]++] = w;
				tandbelow[2] &= sn->all_previous_cells; break;
			}			
			nbelow++;
			goto next;
		}
	}
	sn->possible_cells = ua_ret7p & s->active_cells;
	if (sn->possible_cells) s++;  // switch to next spot
	goto next;
}

void G17B::GoExpand_7_10() {
	SPB03 * sn= &spb_0_15[7];
	t54b12.ntm = 0;
	if (!Expand_7_10()) return;
	p_cpt2g[6]++;
	int nbelow = ntbelow[4] + ntbelow[3] + ntbelow[2] + ntbelow[1] + ntbelow[0];
	p_cpt2g[51] += ntbelow[5];
	p_cpt2g[52] += nbelow;
	if (ntbelow[5] <20)p_cpt2g[53]++;
	//memset(ntbelow, 0, sizeof ntbelow);//7 8 9 10 11 full
	guah54_2.Build2(sn->all_previous_cells, sn->active_cells);
	for (int ib3 = 0; ib3 < genb12.nband3; ib3++)
		genb12.bands3[ib3].guam2done = 0;
	if (nbelow) {
		cout << "Expand_7_10()  [3] " << p_cpt2g[3]
			<< " [4] " << p_cpt2g[4] << " ";
		DumpPotential(0);
	}
	// process 10 clues 
	if (ntbelow[5] ==1) {// get active g2 g3 from guah54_2 direct
		//guah54_2.Dumpall2();
		p_cpt2g[7]++;
		cb3.ncl = 10;
		myb12 = cb3.bf12 = tfull[0];
		cb3.bfcom = sn->all_previous_cells;
		cb3.cbs = tfullcbs[0];
		cb3.g2t = guah54_2.GetG2(cb3.bf12);
		cb3.g3t = guah54_2.GetG3(cb3.bf12);
		//cb3.Dump();
		clean_valid_done = 0;
		for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
			STD_B3& b3 = genb12.bands3[ib3];
			b3.Go(cb3);
			if(clean_valid_done == 2)break;
		}

	}
	else if (ntbelow[5] < 20) {// get active g2 g3 from guah54_3 using common cells 
	}
	else {// work on vectors 
		for (uint32_t iv = 0; iv < ntbelow[5]; iv ++ ) {
			cb3.ncl = 10;
			myb12 = cb3.bf12 = tfull[iv];
			if (t54b12.NotValid(myb12)) continue;
			cb3.bfcom = sn->all_previous_cells;
			cb3.cbs = tfullcbs[iv];
			cb3.g2t = guah54_2.GetG2(cb3.bf12);
			cb3.g3t = guah54_2.GetG3(cb3.bf12);
			//cb3.Dump();
			p_cpt2g[7]++;
			clean_valid_done = 0;
			for (int ib3 = 0; ib3 < genb12.nband3; ib3++) {
				STD_B3& b3 = genb12.bands3[ib3];
				b3.Go(cb3);
				if (clean_valid_done == 2)break;
			}

			//if (iv > 50)break;
		}

	}
}



void G17B::Expand_7_11() {
	SPB03* sl = &spb_0_15[8], * s = sl, * sn;
	T54B12::TUVECT& tuv128 = t54b12.tc128[0];
	uint64_t* twu = tuv128.t;
	if (t54b12.nc128 < 128)t54b12.nc128 = 128;// force adds outside
	*s = spb_0_15[7];	// duplicate 6 for new vector
	s->possible_cells = twu[0];
	s->v = tuv128.v0;// initial nothing done
	//_______ start search clues 7_11
next:	// catch and apply cell in bitfields
	if (GetNextCell(s))	if (--s >= sl)goto next;	else return;
	sn = s + 1;
#ifdef HAVEKNOWN
	if (!((~pk54) & sn->all_previous_cells)) {
		cout << Char54out(sn->all_previous_cells) << " on the path 7_11 ncl=" << sn->ncl << endl;
	}
#endif
	ua_ret7p = 0;
	if (sn->ncl == 11) {// 11 cells
		if (sn->cbs.IsFilt11_17()) goto next;
		GetNextUaAdd(sn);// check adds
		if (ua_ret7p)goto next;
		clean_valid_done = 0;
		p_cpt2g[5]++;
		t54g2.NewVcl(sn->ncl, tclues[s->ncl]);
		t54g3.NewVcl(sn->ncl, tclues[s->ncl]);
		GoAfterExpand(sn);
		goto next;
	}
	if (sn->ncl == 10) {// last step
		if (GetLastAndUa(sn)) // last not empty
			if (!ua_ret7p) goto next;
	}
	else {
		if (sn->v.isNotEmpty())GetNextUa(sn);// first 128
		else GetNextUaAdd(sn);
	}
	if (!ua_ret7p) {// add ua or go below
		if (IsValid7p(sn)) {// got uas to use
			if (sn->ncl == 10)ua_ret7p = anduab12;
		}
		else {// valid below
			clean_valid_done = 1;
			t54g2.NewVcl(sn->ncl, tclues[s->ncl]);
			t54g3.NewVcl(sn->ncl, tclues[s->ncl]);
			GoAfterExpand(sn);//try direct
			GoBelow(sn);
			goto next;
		}
	}
	sn->possible_cells = ua_ret7p & s->active_cells;
	if (sn->possible_cells) { // switch to next spot
		t54g2.NewVcl(sn->ncl, tclues[s->ncl]);
		t54g3.NewVcl(sn->ncl, tclues[s->ncl]);
		s++;
	}
	goto next;
}

void G17B::GoBelow(SPB03* sn) {

}


//____________processing band 3 to the end


#define stack1_54 07007007007007007
int tstack27[3] = { 07007007 ,070070070,0700700700 };
struct CRITB3 {
	uint32_t minix[4],// triplet bf1 bf2 bf3  
		critbf, pairs27, mincount,
		critstack,
		assigned, active,
		ncl, nb3, nmiss;
	BF64 stackf, stackb12, stackmin;
	inline void Init(int ncb12,  CBS& cbsx) {
		memset(this, 0, sizeof(*this));
		ncl = ncb12, nb3 = 17 - ncl;
		memcpy(&stackb12, cbsx.s, sizeof stackb12);
		//minb3 = minb3e;
	}
	inline int  CanNotAddi27(int i27) {
		int bit27 = 1 << i27;
		if (assigned & bit27)return 1;// critical stack common clue assigned
		int imini = i27 / 3, bitmini = 1 << imini, mask = 7 << (3 * imini);
		if (minix[0] & bitmini) {// was a triplet
			minix[0] ^= bitmini;
			minix[1] |= bitmini;
			pairs27 |= bit27;
			critbf ^= bit27;
			return 0;
		}
		if (minix[1] & bitmini) {// was a one pair
			minix[1] ^= bitmini;
			minix[2] |= bitmini;
			pairs27 |= bit27;
			critbf |= mask;
			return 0;
		}

		if (minix[3] & bitmini) {// can not be
			cout << "critb3 can not add to 3 pairs bug" << endl;
			g17b.aigstop = 1;
			return 1;
		}
		// now true add 0->1 or 2_>3
		int  stack = C_stack[i27], bitstack = 1 << stack;
		if (critstack & bitstack)return 1;// not possible

		if (minix[2] & bitmini) {// was 2 pairs
			minix[2] ^= bitmini;
			minix[3] |= bitmini;
			pairs27 |= bit27;
		}
		else {// first pair in the mini row no triplet
			minix[1] |= bitmini;
			pairs27 |= bit27;
			critbf |= (mask ^ bit27);
		}
		stackf.bf.u16[stack]++;
		nmiss--;
		if (!nmiss) {
			critstack = 7;
			active &= critbf;// no more outfield
		}
		else if (stackf.bf.u16[stack] == nb3) {
			critstack |= bitstack;
			active &= (~tstack27[stack]) | critbf;
		}
		return 0;
	}
	inline int  CanNotAddi9(int i9) {
		//can not have a pair active
		int  stack = i9 % 3, bitstack = 1 << stack;
		if (critstack & bitstack)return 1;// not possible
		int mask = 7 << (3 * i9);
		critbf |= mask;
		minix[0] |= 1 << i9;
		stackf.bf.u16[stack]++;
		nmiss--;
		if (!nmiss) {
			critstack = 7;
			active &= critbf;// no more outfield
		}
		else if (stackf.bf.u16[stack] == nb3) {
			critstack |= bitstack;
			active &= (~tstack27[stack]) | critbf;
		}
		return 0;
	}

	inline int Addone(uint32_t i27) {// back 1 if not possible
		int bit27 = 1 << i27, stack = C_stack[i27];
		int imini = i27 / 3, bitmini = 1 << imini, bitstack = 1 << stack;
		assigned |= bit27;
		if (!(bit27 & critbf)) {// clue added outfield
			if (!nmiss)return 1;// not possible
			//if (critstack & bitstack)return 1;// not possible
			nmiss--;
			if (!nmiss) {
				critstack = 7;
				active &= critbf;// no more outfield
			}
			else {
				stackf.bf.u16[stack]++;
				if (stackf.bf.u16[stack] == nb3)critstack |= bitstack;
			}
			return 0;;
		}
		// now add in field within mincount
		if (minix[3] & bitmini) {// 2 clues expected
			critbf ^= bit27; minix[3] ^= bitmini;
			if (!nmiss)active &= critbf;
			return 0;
		}
		if ((minix[2] & bitmini) && (pairs27 & bit27)) {
			// 2 clues if not common clue one more clue
			if (critstack & bitstack)return 1;// not possible
			critbf ^= bit27;	minix[2] ^= bitmini;
			nmiss--;
			if (!nmiss) {
				critstack = 7;
				active &= critbf;// no more outfield
			}
			else {
				stackf.bf.u16[stack]++;
				if (stackf.bf.u16[stack] == nb3)critstack |= bitstack;
			}
			return 0;
		}
		register int mask = ~(7 << (3 * imini));// clear minirow
		critbf &= mask;
		if (!nmiss)active &= critbf;
		if (minix[1] & bitmini) minix[1] ^= bitmini; // one clue expected 		
		else if (minix[2] & bitmini)minix[2] ^= bitmini;//common cell in 2 pairs 
		else if (minix[0] & bitmini)minix[0] ^= bitmini;// triplet
		return 0;
	}
	inline int AddAssign(uint32_t bfa) {// back 1 if not possible
		if (assigned & bfa)return 1; //should never be
		active &= ~bfa; // minimum is to kill new assign
		register int i27, X = bfa;
		while (bitscanforward(i27, X)) {
			X ^= 1 << i27;// clear bit
			if (Addone(i27))return 1;
		}
		return 0;
	}

	void AssignBf2(int bf) {// all or part of the bf2
		minix[2] &= ~bf;// clean the bf2 
		int imini;
		while (bitscanforward(imini, bf)) {
			int  mask = 7 << (3 * imini), m27 = pairs27 & mask;
			bf ^= 1 << imini;// clean the bf
			critbf &= ~mask;// clean the mini row as crit field
			active &= ~mask;// no more clue in this minirow
			pairs27 ^= m27;// clean the 27 pairs to assign
			assigned |= mask ^ m27;// assign the third cell of the mini row
		}
	}
	inline void AssignCritical() {
		active = critbf;
		if (critstack == 7) {
			if (minix[2])  AssignBf2(minix[2]);
			return;
		}
		else if (critstack) {
			for (int i = 0, bit = 1, mask = 07007007; i < 3; i++, bit <<= 1, mask <<= 3)
				if (!(critstack & bit)) active |= mask;
		}
		else active = BIT_SET_27;
		if ((!critstack) || (!minix[2])) return;
		for (int i = 0, bit = 1, mask = 0111; i < 3; i++, bit <<= 1, mask <<= 1) {
			if (critstack & bit) {
				active &= ~(stack1_54 << 3 * i); //  (stack154 << 3*i)active must be in critbf for this; stack
				if (minix[2] & mask)	AssignBf2(minix[2] & mask);

			}
		}
		active |= critbf;// be sure to keep it

	}

	inline void MinVmini(int imini, int vmini) {
		uint32_t shift = 3 * imini, bit9 = 1 << imini, cmin = 1,
			mask = 7 << shift, vminishift = vmini << shift;
		if (vmini == 8) {
			minix[0] |= bit9;// mini triplet
			critbf |= mask;
		}
		else {
			pairs27 |= vminishift;
			uint32_t cc = _popcnt32(vmini);
			minix[cc] |= bit9;
			if (cc > 1)critbf |= mask;
			else critbf |= (mask ^ vminishift);
			if (cc == 3) 	cmin = 2;
		}
		mincount += cmin;
		stackmin.bf.u16[imini % 3] += cmin;
	}
	inline int Stackx() {
		register uint32_t x = nb3;
		if (stackf.bf.u16[0] > x)return 1;
		if (stackf.bf.u16[1] > x)return 1;
		if (stackf.bf.u16[2] > x)return 1;
		if (stackf.bf.u16[0] == x)critstack |= 1;
		if (stackf.bf.u16[1] == x)critstack |= 2;
		if (stackf.bf.u16[2] == x)critstack |= 4;
		return 0;
	}
	inline void SetStackf(uint64_t bf) {
		stackf.bf.u64 = stackmin.bf.u64 + stackb12.bf.u64;
		nmiss = nb3 - mincount;
		if (!nmiss) critstack = 7;
	}
	inline int GetToAss() { return (nb3 - _popcnt32(assigned)); }
	void Status(const char* lib) {
		cout << lib << "critical Status mincount =" << mincount << " nmiss=" << nmiss
			<< " critstack=" << critstack << endl;
		cout << Char27out(critbf) << " critical bf" << endl;
		cout << Char27out(pairs27) << " pairs 27" << endl;
		cout << Char27out(assigned) << " assigned" << endl;
		cout << Char27out(active) << " active" << endl;
		if (minix[1])cout << Char9out(minix[1]) << "     minis bf1" << endl;
		if (minix[2])cout << Char9out(minix[2]) << "     minis bf2" << endl;
		if (minix[3])cout << Char9out(minix[3]) << "     minis bf3" << endl;
		if (minix[0])cout << Char9out(minix[0]) << " mini triplets" << endl << endl;
		cout << "b12 stacks " << stackb12.bf.u16[0] << stackb12.bf.u16[1] << stackb12.bf.u16[2] << endl;
		cout << "minis stacks " << stackmin.bf.u16[0] << stackmin.bf.u16[1] << stackmin.bf.u16[2] << endl;
		cout << "final stacks " << stackf.bf.u16[0] << stackf.bf.u16[1] << stackf.bf.u16[2] << endl;
	}

}scritb3;

void G17B::GoAfterExpand(SPB03* sn, uint32_t nadd) {
	if (aigstop) return;
#ifdef HAVEKNOWN
	if (pk54 == sn->all_previous_cells) {
		cout << Char54out(sn->all_previous_cells) << " expected GoAfterExpand nadd=" << endl;
	}
#endif
	p_cpt2g[6]++;

	ntuaddb12 = nadd;// used only if add mode
	nclgo = sn->ncl;
	ncluesb3 = 18 - nclgo;
	myb12f = sn->all_previous_cells;

	/*



	scritb3.SetStackf(Bf);// stackf; nmisss; critstack if nmiss=0
	if (scritb3.Stackx())return;
	scritb3.AssignCritical();
	ntw3 = 0; // temporary before  split in/out
	g_256h.ApplyMore(&tclues[6], nclgo - 6,
		tw3, ntw3, scritb3.assigned);
	register uint32_t  critbf = scritb3.critbf,
		nmiss = scritb3.nmiss, n1 = 0, n2 = 0;
	if (!nmiss) {// first out is "dead"
		for (uint32_t i = 0; i < ntw3; i++)
			if (!(tw3[i] & critbf)) return;
	}
	else if (nmiss == 1) {//use "and" as out
		register uint32_t andout = scritb3.active & (~critbf);
		for (uint32_t i = 0; i < ntw3; i++) {
			register uint32_t U = tw3[i];
			if (U & critbf) tw3[n1++] = U;
			else { andout &= U; n2 = 1; }
		}
		if (n2) {
			if (!andout)return;// no possibility to add
			tw3[n1++] = andout; // dummy ua hitting all "out"
		}
		ntw3 = n1;
	}
	ntaddgob3 = 0;
#ifdef HAVEKNOWN
	if (pk54 == sn->all_previous_cells) cout << "call buildb2" << endl;
#endif
	if (BuildB2Table()) return;
	int nass = _popcnt32(scritb3.assigned);
	if (nass == scritb3.nb3) {// already ok to go
		if (ntw3_2)return;
		if (scritb3.BadPerm()) return;
		if (!clean_valid_done) {
			clean_valid_done = 1;
			if (IsValidB12(nclgo)) { clean_valid_done = 2; return; }
		}
		if (!IsValidB3(scritb3.assigned)) Out17(scritb3.assigned);
		return;
	}
#ifdef HAVEKNOWN
	if (pk54 == myb12f) cout << "call BuildFinalTable" << endl;
#endif
	BuildFinalTable();
#ifdef HAVEKNOWN
	//if (pk54 == sn->all_previous_cells) aigstop = 1;
#endif
}*/

}

int G17B::IsValid_myb12() {
	nclues = 0;
	int cell;
	register uint64_t U = myb12;
	//cout << Char54out(U) << " isvalid [7] " << p_cpt2g[7] << endl;
	while (bitscanforward64(cell, U)) {
		U ^= (uint64_t)1 << cell;
		tclues[nclues++] = cell;
	}
	if (zh2b[1].IsValid(tclues, nclues)) {
		anduab12 = ~0;
		register uint64_t ua54;
		for (uint32_t i = 0; i < zh2gxn.nua; i++) {
			register uint64_t ua = zh2gxn.tua[i], cc = _popcnt64(ua);
			ua54 = (ua & BIT_SET_27) | ((ua & BIT_SET_B2) >> 5);
			anduab12 &= ua54;
			if (0) {
				int digs = 0;
				for (uint64_t i = 0, bit = 1; i < 54; i++, bit <<= 1)
					if (ua54 & bit) digs |= 1 << grid0[i];
				cout << Char54out(ua54) << " to add 12 size " << cc;
				cout << " "<<Char9out(digs)<<" [7]"<< p_cpt2g[7] << endl;
			}
			if (cc < 22) {
				t54b12.AddA(ua54);
				t54b12.AddB(ua54);
				t54b12.AddC(ua54);
				t54b12.AddM(ua54);
			}
		}
		ua_ret7p = ua54;// return last (smaller)
		return 1;
	}
	return 0;
}
uint32_t G17B::IsValidB3(uint32_t bf) {

	if (zhou[0].CallCheckB3(tclues, nclues, bf)) {
		anduab3 = BIT_SET_27;
		stopexpandb3 = 0;
		//cout << "\t\t\t\t   add b3 nadd="<< zhgxn.nua << endl;
		for (uint32_t iadd = 0; iadd < zhgxn.nua; iadd++) {
			BF128 w = zhgxn.tua[iadd];
			int cc = _popcnt32(w.bf.u32[2]);
			if (!cc) {
				cout << " bug no b3 IsValidB3 [7]" << p_cpt2g[7]
					<< "   [8]" << p_cpt2g[8] << endl;
				aigstop = 1;	return 1;
			}
			register uint32_t ua = w.bf.u32[2];
			t3_2[nt3_2++] = ua;
			anduab3 &= ua;
			register uint64_t U = w.bf.u64[0];
			U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54


			if (cc <= 6 && cc>1) {
				if (cc > 3)myband3->Addguam(w);
				else {
					//cout << Char2Xout(w.bf.u64[0]) << " | ";
					//cout << Char27out(w.bf.u32[2]) << " add b3" << endl;
					if (cc == 2) {
						int i81 = myband3->GetI81_2(w.bf.u32[2]);
						guah54.Add2(U, i81); guah54_2.Add2(U, i81);
						cb3.g2t.setBit(i81);
						p_cpt2g[20]++;
					}
					else {
						int i81 = myband3->GetI81_3(w.bf.u32[2]);
						guah54.Add3(U, i81); guah54_2.Add2(U, i81);
						cb3.g3t.setBit(i81);
						p_cpt2g[21]++;
					}
				}

			}
			/*
					else if(cc3==2  &&  ccua>10) { // must build g2
					int i81=myband3->GetI81_2(w.bf.u32[2]);
					guah.tg2[i81].AddIf(wd.bf.u64[0]);
				}
				else if (cc3 == 3  &&  ccua>9) { // must build g3
					int i81 = mb3.GetI81_3(wd.bf.u32[2]);
					guah.tg3[i81].AddIf(wd.bf.u64[0]);
				}
		uint64_t cc0 = _popcnt64(w.bf.u64[0]);
			taddgob3[ntaddgob3++] = ua;// more for a given b1+2
			BuildGua(w, cc);
			g_256h.Add128(w, cc, 1);//clean level new if possible
			if (cc < 4) {
				if ((cc == 2 && scritb3.CanNotAddi27(w.bf.u32[2])) ||
					(cc == 3 && scritb3.CanNotAddi9(w.bf.u32[2]))) {
					stopexpandb3 = 1;
				}
			}
			if (cc0 > 20 || cc > 6)continue;
			chunkh.Add128(w, cc);
			if (cc < 3)p_cpt2g[35]++; else p_cpt2g[36]++;
			*/
		}
		return 1;
	}
	else return 0;
}

void STD_B3::Go(CALLBAND3& cb3) {
	p_cpt2g[8]++;
	//g.gsocket2.Print("g2"); g.gsocket3.Print("g3");
	BF128 g2 = cb3.g2t & g.gsocket2,	g3 = cb3.g3t & g.gsocket3;
	//g2.Print("g2&"); g3.Print("g3&");
	register uint32_t  vmini, Mg2 = 0, Mg3 = 0;
	int x;
	while ((x = g2.getFirst128()) >= 0) {
		//cout << Char27out(g.ua2bit27[x]) << " x=" << x << endl;
		g2.clearBit(x);		Mg2 |= g.ua2bit27[x];	}
	while ((x = g3.getFirst128()) >= 0) {
		//cout << Char9out(g.ua3_imini[x]) << " x=" << x << endl;
		g3.clearBit(x);		Mg3 |= 1<<g.ua3_imini[x];	}
	//cout << Char27out(Mg2) << " Mg2" << endl;
	//cout << Char9out(Mg3) << " Mg3" << endl;
	scritb3.Init(cb3.ncl, cb3.cbs);
	for (int imini = 0; imini < 9; imini++, Mg2 >>= 3, Mg3 >>= 1) {
		if (!(vmini = Mg2 & 7))	if (Mg3 & 1)vmini = 8;
		if (vmini)scritb3.MinVmini(imini, vmini);
	}
	if ((int)scritb3.mincount >   scritb3.nb3) return;
	//cout << band << " b3 to process" << endl;
	//scritb3.Status(" aaaa  ");
	scritb3.SetStackf(cb3.bf12);// stackf; nmisss; critstack if nmiss=0
	if (scritb3.Stackx())return;// never > nb3 in stack
	scritb3.AssignCritical();// can assign 2 pairs in critical stacks
	//scritb3.Status(" bbbbbb  ");
	if (!guam2done) BuildGuam2(cb3.bfcom);

	// now get still valid uas band 3 not gua2s gua3s
	g17b.nt3 = 0;
	register uint64_t F = cb3.bf12;
	register uint32_t As = scritb3.assigned,Ac= scritb3.active;
	if (0) {
		cout << Char54out(F) << " F" << endl;
		cout << Char27out(As) << " As" << endl;
		cout << Char27out(Ac) << " Ac" << endl;
	}
	{
		register uint32_t V = scritb3.pairs27,i27;
		while (bitscanforward64(i27, V)) {
			V ^= 1 << i27;
			g17b.t3[g17b.nt3++] = g.pat2[i_27_to_81[i27]];
		}
		V = scritb3.minix[0];// triplets
		if (V) {// add guas3 if any (usually 0)
			for (int i = 0; i < 9; i++)if (V & (1 << i))
				g17b.t3[g17b.nt3++] = 7 << (3 * i);
		}

	}
	for (uint32_t i = 0; i < nua; i++) {// now band 3 uas
		register uint32_t U = tua[i];
		if (U & As) continue;// assigned in bands 3
		//cout << Char27out(U) << " Ub3" << endl;
		if (!(U &= Ac)) return; // dead branch
		g17b.t3[g17b.nt3++] = U;
	}

	for (uint32_t i = 0; i < ntguam2; i++) {
		GUAM gm= tguam2[i];
		if (gm.bf12 & F) continue;// assigned in bands 1+2
		register uint32_t U = gm.bf3;
		if (U & As) continue;// assigned in bands 3
		//cout << Char27out(U) << " Uguam"  << endl;
		if (!(U &= Ac)) return; // dead branch
		g17b.t3[g17b.nt3++] = U;
	}
	g17b.myband3 = this;
	g17b.GoB3CleanOne();



}

void G17B::GoB3CleanOne() {// assign all singles in uas
	while (1) {
		//for (int i = 0; i < (int)nt3; i++) {
			//register uint32_t U = t3[i];
			//cout << Char27out(t3[i]) << " i=" << i << endl;
		//}
		int is1 = 0;
		for (int i = 0; i < (int)nt3; i++) {
			register uint32_t U = t3[i];
			if (_popcnt32(U) == 1)  is1 |= U;			
		}
		if (is1) {
			//cout << Char27out(is1) << " is1" << endl;
			if (scritb3.AddAssign(is1)) return ;// conflict
			register uint32_t AC = scritb3.active,	F = scritb3.assigned;
			int n = nt3; nt3 = 0;
			for (int i = 0; i < n; i++) {
				register uint32_t U = t3[i];
				if (!(U & F)) {
					if (!(U &= AC)) return; // dead branch
					t3[nt3++] = U ;
				}
			}
		}
		else break;
	}
	memcpy(&grid0[54], myband3->band0, sizeof genb12.bands3[0].band0);

	if (_popcnt32(scritb3.assigned) == scritb3.nb3) {
		if (!clean_valid_done) {
			clean_valid_done = 1;
			if (IsValid_myb12()) { clean_valid_done = 2; return; }
		}
		if (clean_valid_done == 2) return;
		//cout << Char54out(myb12)<< " [7]"<< p_cpt2g[7]  << endl;
		//cout << myband3->band << " filled cleanone " << endl;
		//cout << Char27out(scritb3.assigned) << " to check valid" << endl;
		//cout << " is valid  b12" << endl;
		if (IsValidB3(scritb3.assigned)) return;
		cout << " is valid 17" << endl;
		//_________________________   valid 18 
		Out17(scritb3.assigned);
		return;
	}
	switch (scritb3.nmiss) {
	case 0: GoB3Miss0(); return;
	case 1: GoB3Miss1(); return;
	default:GoB3MissMore(); return;
	}
}
void G17B::GoB3Miss0() {// no out field here assigned expected
	// store final  and check redundancy
	p_cpt2g[10]++;
	nt3_2 = 0;
	for (uint32_t i = 0; i < nt3; i++) {
		register uint32_t U = t3[i], nu = ~U;
		for (uint32_t j = 0; j < nt3_2; j++) {
			if (!(t3_2[j] & nu)) { U = 0; break; }
		}
		if (U)t3_2[nt3_2++] = U;
	}
	if (0) {
		cout << Char54out(myb12) << endl;
		cout << myband3->band << " b3 after cleanone nt3=" << nt3 << endl;
		scritb3.Status(" after clean one ");
		cout << "final table" << endl;
		for (int i = 0; i < (int)nt3_2; i++) {
			register uint32_t U = t3_2[i];
			cout << Char27out(U) << " i=" << i << endl;
		}
	}
	int ntoass = scritb3.GetToAss();
	if (ntoass > 1) {
		ExpandB3Direct(ntoass);
		return;
	}
	// only one to ass must hit all uas
	register int A= t3_2[0];
	for (int i = 1; i < (int)nt3_2; i++) 
		if (!(A &= t3_2[i])) return;
	if (!clean_valid_done) {
		clean_valid_done = 1;
		if (IsValid_myb12()) { clean_valid_done = 2; return; }
	}
	register int cell;
	while (bitscanforward(cell, A)) {
			CRITB3 critb3 = scritb3;
			register uint32_t bit = 1 << cell;
			A ^= bit; //clear bit
			if (critb3.Addone(cell)) 	continue;
			if (IsValidB3(critb3.assigned)) {
				if (stopexpandb3) return;
				A &= anduab3;
			}
			else 	Out17(critb3.assigned);
	}
	

}
void G17B::GoB3Miss1() {// if outfield must have one out 
	p_cpt2g[11]++;

	int n2 = 0;
	{// store final  get outfield and check redundancy
		nt3_2 = 0;
		register uint32_t  critbf = scritb3.critbf,
			nmiss = scritb3.nmiss, n1 = 0;
		register uint32_t andout = scritb3.active & (~critbf);
		for (uint32_t i = 0; i < nt3; i++) {
			register uint32_t U = t3[i], nu = ~U;
			if (!(U & critbf)) { andout &= U; n2 = 1; continue; }
			for (uint32_t j = 0; j < nt3_2; j++) {
				if (!(t3_2[j] & nu)) { U = 0; break; }
			}
			if (U)t3_2[nt3_2++] = U;
		}
		if (n2) {
			if (!andout)return;// no possibility to add
			t3_2[nt3_2++] = andout; // dummy ua hitting all "out"
		}
	}
	if (!clean_valid_done) {// stop if not a valid band 1+2
		clean_valid_done = 1;
		if (IsValid_myb12()) { clean_valid_done = 2; return; }
	}
	if ((!n2) || (!scritb3.minix[2])) {// nothing to assign
		ExpandB3Direct(scritb3.GetToAss());
		return;
	}

	if (0) {
		cout << Char54out(myb12) << " [7]" << p_cpt2g[7] << " cdone=" << clean_valid_done << endl;
		cout << myband3->band << " b3 after cleanone nt3=" << nt3 << endl;
		scritb3.Status(" after clean one ");
		cout << "final table n2(outfiel)=" << n2 << endl;
		for (int i = 0; i < (int)nt3_2; i++) {
			register uint32_t U = t3_2[i];
			cout << Char27out(U) << " i=" << i << endl;
		}
	}


	// now all critical can assign bf2
	scritb3.critstack = 7;
	scritb3.AssignCritical();
	scritb3.active |= t3_2[nt3_2 - 1];// outfield lost to reenter 
	//scritb3.Status("after assign critical");
	// reduce uas and look for fresh singles
	while (1) {
		int is1 = 0;
		register uint32_t AC = scritb3.active, F = scritb3.assigned;
		int n = nt3_2; nt3_2 = 0;
		for (int i = 0; i < n; i++) {
			register uint32_t U = t3_2[i];
			if (!(U & F)) {
				if (!(U &= AC)) return; // dead branch
				if (_popcnt32(U) == 1)  is1 |= U;
				//cout << Char27out(U) << " loop nt3_2=" << nt3_2 << endl;
				t3_2[nt3_2++] = U;
			}
		}

		if (is1) {
			//cout << Char27out(is1) << " is1" << endl;
			if (scritb3.AddAssign(is1)) return;// conflict
		}
		else break;
	}
	//scritb3.Status("after assign singles");
	ExpandB3Direct(scritb3.GetToAss());

}
void G17B::GoB3MissMore() {// no out field analysis here 
	// store final  and check redundancy
	p_cpt2g[12]++;
	nt3_2 = 0;
	for (uint32_t i = 0; i < nt3; i++) {
		register uint32_t U = t3[i], nu = ~U;
		for (uint32_t j = 0; j < nt3_2; j++) {
			if (!(t3_2[j] & nu)) { U = 0; break; }
		}
		if (U)t3_2[nt3_2++] = U;
	}
	if (nt3_2 < 25) {
		ExpandB3Direct(scritb3.GetToAss());
		return;
	}
	p_cpt2g[14]++;


	if (p_cpt2g[14]<10) {
		cout << Char54out(myb12) <<" [7] "<< p_cpt2g[7] << endl;
		cout << myband3->band << " b3 after cleanone nt3=" << nt3 << endl;
		scritb3.Status(" after clean one ");
		cout << "final table" << endl;
		for (int i = 0; i < (int)nt3_2; i++) {
			register uint32_t U = t3_2[i];
			cout << Char27out(U) << " i=" << i << endl;
		}
	}
}

/*
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
	}*/

void G17B::GoB3small(STD_B3* mb3) {// most cases 
	// sort and pack
	BF128 vsize[26];// vsize [25] is flag for out field
	memset(vsize, 0, sizeof vsize);
	register uint32_t InField = scritb3.critbf, Andout = ~0;
	for (int i = 0; i < (int)nt3; i++) {
		register uint32_t U = t3[i];
		//cout << Char27out(t3[i]) << "i=" << i << endl;
		int cc = _popcnt32(U);
		if (!(U & InField)) {
			vsize[25].setBit(i);
			Andout &= U;
		}
		if (cc > 24)cc = 24;//should never be
		vsize[cc].setBit(i);
	}


}

/*
ntuaddb12 = nadd;// used only if add mode


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

	ntw3 = 0; // temporary before  split in/out
	g_256h.ApplyMore(&tclues[6], nclgo - 6,
		tw3, ntw3, scritb3.assigned);
	register uint32_t  critbf = scritb3.critbf,
		nmiss = scritb3.nmiss, n1 = 0, n2 = 0;
	if (nmiss == 1) {//use "and" as out
		register uint32_t andout = scritb3.active & (~critbf);
		for (uint32_t i = 0; i < ntw3; i++) {
			register uint32_t U = tw3[i];
			if (U & critbf) tw3[n1++] = U;
			else { andout &= U; n2 = 1; }
		}
		if (n2) {
			if (!andout)return;// no possibility to add
			tw3[n1++] = andout; // dummy ua hitting all "out"
		}
		ntw3 = n1;
	}



*/
void G17B::GoB3Big(STD_B3* mb3) {// most cases 
	// sort and pack
	BF128 vsize[26];// vsize [25] is flag for out field
	memset(vsize, 0, sizeof vsize);
	register uint32_t InField = scritb3.critbf, Andout = ~0;
	if (1) {
		for (int i = 0; i < (int)nt3; i++) {
			register uint32_t U = t3[i];
			cout << Char27out(t3[i]) << "i=" << i << endl;

		}
	}
	cout << "need more than 128 uas" << endl;

}
void G17B::ExpandB3Direct(int ntoass) {
	uint64_t limspot = ntoass - 1, limm1 = limspot - 1;
	//cout << "b3direct lims " << limspot << " " << limm1<<" clean_valid_done="<< clean_valid_done << endl;
	struct SPB {
		CRITB3 critb3;
		uint32_t  possible_cells, indtw3;
	}spb[12];

	register SPB* s, * sn;
	register uint64_t ispot;
	s = spb;
	memset(s, 0, sizeof spb[0]);
	s->critb3 = scritb3;
	s->indtw3 = 0;// initial first ua
	s->possible_cells = t3_2[0];
	//cout << Char27out(s->possible_cells) << " start" << endl;
next:
	ispot = s - spb;
	{// catch and apply cell in bitfields
		register int cell;
		register uint32_t p = s->possible_cells;
		if (!p)if (--s >= spb)goto next; else return;
		bitscanforward(cell, p);
		register uint32_t bit = 1 << cell;
		s->possible_cells ^= bit;
		s->critb3.active ^= bit;
		sn = s + 1; *sn = *s;
		if (sn->critb3.Addone(cell)) 	goto next;
	}
	if (ispot < limm1) {// here max 16 clues never valid b3
		//if (locdiag) cout << " not last step s->indtw3="<< s->indtw3 << endl;
		register uint32_t F = sn->critb3.assigned;
		for (uint32_t i = s->indtw3 + 1; i < nt3_2; i++) {
			register uint32_t U = t3_2[i];
			if (!(U & F)) {
				U &= s->critb3.active;
				if (!U)goto next;//dead branch
				sn->possible_cells = U;
				//cout << Char27out(U) << " next ispot="<<ispot << endl;
				sn->indtw3 = i;
				s++; // switch to next spot
				goto next;
			}
		}
		//cout << "check before limit or n to ass=1 and one ua" << endl;
		// no ua available must check( in 18 mode can not be valid)
		if (!clean_valid_done) {
			clean_valid_done = 1;
			if (IsValid_myb12()) { clean_valid_done = 2; return; }
		}
		if (IsValidB3(sn->critb3.assigned)) {
			if (stopexpandb3) return;
			s->possible_cells &= anduab3;// same spot hitting  new ua
			sn->possible_cells = anduab3;
			s = sn;
			goto next;
		}
		cout << "bug seen valid 16 or less [5] " << endl;
		aigstop = 1;
		return;
	}

	// this is the last step must hit all pending uas
	{ // find next cells hitting all uas
		int aig = 1;
		register uint32_t andw = sn->critb3.active;
		register uint32_t F = sn->critb3.assigned;
		//cout <<Char27out(andw) << " active entry last" << endl;
		//cout << Char27out(F) << " assigned entry last" << endl;
		{// and still valid   uas
			for (uint32_t i = s->indtw3 + 1; i < nt3_2; i++) {
				register uint32_t U = t3_2[i];
				if (!(U & F)) { // not hit ua
					if (!(andw &= U))goto next;//dead branch
					aig = 0;
				}
			}
		}
		// no more ua or "and all uas" not empty
		if (!clean_valid_done) {
			clean_valid_done = 1;
			if (IsValid_myb12()) { clean_valid_done = 2; return; }
		}
		if (aig) {// no ua could be a 17 valid in 18 mode 
			if (IsValidB3(F)) {
				if (stopexpandb3) return;
				andw &= anduab3;
			}
			else { Out17(F);	goto next; }// this is a valid 17
		}
		//cout << Char27out(andw) << " last ispot=" << ispot << endl;

		register int cell;
		while (bitscanforward(cell, andw)) {
			CRITB3 critb3 = sn->critb3;
			register uint32_t bit = 1 << cell;
			andw ^= bit; //clear bit
			if (critb3.Addone(cell)) 	continue;
			if (IsValidB3(F | bit)) {
				if (stopexpandb3) return;
				andw &= anduab3;
			}
			else 	Out17(F | bit);
		}
		goto next;
	}
	goto next;// never here
}


/*


void GCHK::Expand_7_12() {
	SPB03* sl = &spb_0_15[8], * s = sl, * sn;
	if (aigstop) return;
	TUASB12::TUAVECT& tuv128 = tuasb12.tv128h.tvc128[0];
	uint64_t* twu = tuv128.t;
	if (tuasb12.tv128h.nc128 < 128)tuasb12.tv128h.nc128 = 128;// force adds outside
	*s = spb_0_15[7];	// duplicate 6 for new vector
	s->possible_cells = twu[0];
	s->v = tuv128.v0;// initial nothing done
	Init7p_guas();// do init 7 clues

	//_______ start search clues 7_12
next:	// catch and apply cell in bitfields
	if (GetNextCell(s))	if (--s >= sl)goto next; else return;
	sn = s + 1;
#ifdef HAVEKNOWN
	if (!((~pk54) & sn->all_previous_cells)) {
		cout << Char54out(sn->all_previous_cells) << " on the path 7_12 ncl=" << sn->ncl << endl;
	}
#endif
	ua_ret7p = 0;
	if (sn->ncl == 12) {// 12 cells
		if (sn->cbs.IsFilt12()) goto next;
		p_cpt2g[5]++;
		GetNextUaAdd(sn);// check adds
		if (ua_ret7p)goto next;
		g_256h.NewVcl(sn->ncl, tclues[s->ncl]);
		clean_valid_done = 0;
		GoAfterExpand(sn);
		if (aigstop) return;
		goto next;
	}
	ua_ret7p = 0;
	if (sn->ncl == 11) {// last step
		if (GetLastAndUa(sn)) {// first 128
			if (!ua_ret7p) goto next;
		}
	}
	else {
		if (sn->v.isNotEmpty())		GetNextUa(sn);// first 128
		else GetNextUaAdd(sn);
	}
	if (!ua_ret7p) {
#ifdef HAVEKNOWN
		if (!((~pk54) & sn->all_previous_cells)) {
			cout << Char54out(sn->all_previous_cells) << " after on the path 7_12 ncl=" << sn->ncl << endl;
		}
#endif
		if (sn->ncl < 11) {
			if (!IsValid7p(sn)) {// valid below
				g_256h.NewVcl(sn->ncl, tclues[s->ncl]);
				clean_valid_done = 1;
				GoAfterExpand(sn);//try direct
				GoBelow(sn);		goto next;
			}
		}
		else {// 11 clues try direct first
			if (sn->cbs.IsFilt11())goto next;
			g_256h.NewVcl(sn->ncl, tclues[s->ncl]);
			clean_valid_done = 0;
			GoAfterExpand(sn);//try direct
			if (clean_valid_done == 2)  // bach not validb12
				ua_ret7p = anduab12;
			else { GoBelow(sn); goto next; }
		}
	}
	sn->possible_cells = ua_ret7p & s->active_cells;
	if (sn->possible_cells) {
		g_256h.NewVcl(sn->ncl, tclues[s->ncl]);
		s++; // switch to next spot
	}
	goto next;
}
*/


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


//=========brute force specific to this
int ZHOU::CallCheckB3(uint32_t* t, int n, uint32_t bf, int nogo) {// 17 search mode
	zh_g.diag = nogo;
	memcpy(this, zhoustart, sizeof zhoustart);
	misc.SetAll_0();
	BF128 dca[9];
	memset(dca, 0, sizeof dca);
	int digitsbf = 0;
	for (int icell = 0; icell < n; icell++) {
		int cell = t[icell], digit = g17b.grid0[cell];
		int xcell = C_To128[cell]; // the cell value in 3x32 of a 128 bits map
		digitsbf |= 1 << digit;
		Assign(digit, cell, xcell);
		dca[digit].Set(xcell);
	}
	{
		uint32_t cc;
		register int x = bf;
		while (bitscanforward(cc, x)) {
			x ^= 1 << cc; //clear bit
			int cell = cc + 54, digit = g17b.grid0[cell];
			digitsbf |= 1 << digit;
			int xcell = cc + 64; // the cell value in 3x32 of a 128 bits map
			Assign(digit, cell, xcell);
			dca[digit].Set(xcell);
		}
	}
	zh_g2.s17_b3_mini = 1;
	BF128 w = cells_unsolved;
	w.bf.u32[3] = ~0;// keep rowunsolved settled
	for (int i = 0; i < 9; i++)  FD[i][0] &= w | dca[i];
	//__________end assign last lot start solver
	zhgxn.nua = 0;
	zh_g.go_back = 0;	zh_g.nsol = 0; zh_g.lim = 1;// modevalid is set to  1
	//ImageCandidats();
	int ir = Full17Update();
	if (ir == 2) return 0;// solved can not be multiple
	Guess17(0);
	return zhgxn.nua;
}

int ZHOU::PartialInitSearch17(uint32_t* t, int n) {
	zh_g2.digitsbf = 0;
	memset(zh_g2.Digit_cell_Assigned, 0, sizeof zh_g2.Digit_cell_Assigned);
	memcpy(this, zhoustart, sizeof zhoustart);
	for (int icell = 0; icell < n; icell++) {
		int cell = t[icell], digit = zh_g2.grid0[cell];
		int xcell = C_To128[cell]; // the cell value in 3x32 of a 128 bits map
		zh_g2.digitsbf |= 1 << digit;
		if (FD[digit][0].Off(xcell))  return 1;// check not valid entry
		Assign(digit, cell, xcell);
		zh_g2.Digit_cell_Assigned[digit].Set(xcell);
	}
	//cout << "ok init" << endl;
	BF128 w = cells_unsolved;
	w.bf.u32[3] = ~0;// keep rowunsolved settled
	for (int i = 0; i < 9; i++)  FD[i][0] &= w | zh_g2.Digit_cell_Assigned[i];
	return 0;
}

int ZHOU::Apply17SingleOrEmptyCellsB12() {
	zh_g.single_applied = 0;
	// here  singles and empty cells till 4 cells searched
	register uint64_t R1, R2, R3, R4;
	{
		register uint64_t* P = FD[0][0].bf.u64, M = *P;
		R1 = M;
		P += 4; M = *P;	                            R2 = R1 & M;  R1 |= M;
		P += 4; M = *P;              R3 = R2 & M;   R2 |= R1 & M; R1 |= M;
		P += 4; M = *P; R4 = R3 & M;  R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 4; M = *P;	R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 4; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 4; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 4; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 4; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
	}
	if (cells_unsolved.bf.u64[0] & (~R1)) 	return 1; // empty cells
	R1 &= ~R2; // now true singles
	R1 &= cells_unsolved.bf.u64[0]; // these are new singles
	if (R1) {// singles to apply
		uint32_t xcell;
		while (bitscanforward64(xcell, R1)) {
			R1 ^= (uint64_t)1 << xcell;
			uint32_t cell = From_128_To_81[xcell];
			for (int idig = 0; idig < 9; idig++) {
				if (FD[idig][0].On(xcell)) {
					Assign(idig, cell, xcell);
					goto nextr1;
				}
			}
			return 1; // conflict with previous assign within this lot
		nextr1:;
		}
		zh_g.single_applied = 1;
		return 0;
	}
	// no single store apply  pair in priority ??
	R2 &= ~R3; // now true singles
	if (!R2) {
		R3 &= ~R4;
		if (R3) R2 = R3;
		else R2 = R4;
	}
	bitscanforward64(zh_g2.xcell_to_guess, R2);
	return 0;
}
int ZHOU::Apply17SingleOrEmptyCellsB3() {
	zh_g.single_applied = 0;
	// here  singles and empty cells till 4 cells searched
	register uint32_t R1, R2, R3, R4;
	{
		register uint32_t* P = &FD[0][0].bf.u32[2], M = *P;
		R1 = M;
		P += 8; M = *P;	                            R2 = R1 & M;  R1 |= M;
		P += 8; M = *P;              R3 = R2 & M;   R2 |= R1 & M; R1 |= M;
		P += 8; M = *P; R4 = R3 & M;  R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 8; M = *P;	R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 8; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 8; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 8; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
		P += 8; M = *P; R4 |= R3 & M; R3 |= R2 & M; R2 |= R1 & M; R1 |= M;
	}
	if (cells_unsolved.bf.u32[2] & (~R1)) 	return 1; // empty cells
	R1 &= ~R2; // now true singles
	R1 &= cells_unsolved.bf.u32[2]; // these are new singles
	if (R1) {// singles to apply in band 3
		uint32_t xcell;
		while (bitscanforward(xcell, R1)) {
			R1 ^= (uint64_t)1 << xcell;
			uint32_t cell = xcell + 54;
			xcell += 64;
			for (int idig = 0; idig < 9; idig++) {
				if (FD[idig][0].On(xcell)) {
					Assign(idig, cell, xcell);
					goto nextr1;
				}
			}
			return 1; // conflict with previous assign within this lot
		nextr1:;
		}
		zh_g.single_applied = 1;
		return 0;
	}
	// no single store apply  pair in priority ??
	R2 &= ~R3; // now true pairs
	if (!R2) {
		R3 &= ~R4; // now true tripletss
		if (R3) R2 = R3;
		else R2 = R4;
	}
	bitscanforward(zh_g2.xcell_to_guess, R2);
	zh_g2.xcell_to_guess += 64;
	return 0;
}
int ZHOU::Full17Update() {
	if (zh_g.go_back) return 0;
	while (1) {
		if (!Update()) return 0; // game locked in update
		//if(diag &&cells_unsolved.bf.u32[2]==0)ImageCandidats();
		if (!Unsolved_Count()) return 2;
		if (cells_unsolved.bf.u32[2]) {// fill B3 first
			if (Apply17SingleOrEmptyCellsB3())	return 0; //  empty cell or conflict singles in cells
		}
		else {
			if ((!ISFALSEON))return 0;
			if (Apply17SingleOrEmptyCellsB12())	return 0;
		}
		if (!zh_g.single_applied)	break;
	}
	return 1;
}
void ZHOU::Compute17Next(int index) {
	int ir = Full17Update();
	//if (zh_g.diag) {
		//cout << "index=" << index << endl;
		//ImageCandidats();
	//}
	if (!ir) return;// locked 
	if (ir == 2) {//solved
		if (index) {// store false as ua
			BF128  wua;
			int* sol = g17b.grid0;
			wua.SetAll_0();
			for (int i = 0; i < 81; i++) {
				int d = sol[i];
				if (FD[d][0].Off_c(i))	wua.Set_c(i);
			}
			if (wua.isNotEmpty()) {
				//if (zh_g.diag) {
					//cout << "zhgxn.nua=" << zhgxn.nua << endl;
					//wua.Print3(" ");
				//}
				int cc = _popcnt32(wua.bf.u32[2]);
				if ((!zhgxn.nua) || cc < 3)
					zhgxn.tua[zhgxn.nua++] = wua;
				if (cc < 3 || zhgxn.nua>5)	zh_g.go_back = 1;
			}
		}
		return;
	}
	Guess17(index);// continue the process
}
void ZHOU::Guess17(int index) {
	if (zh_g.go_back) return;
	int xcell = zh_g2.xcell_to_guess,
		cell = From_128_To_81[xcell],
		digit = zh_g2.grid0[cell];
	// true first if possible
	if (FD[digit][0].On(xcell)) {
		ZHOU* mynext = (this + 1);
		*mynext = *this;
		mynext->SetaCom(digit, cell, xcell);
		mynext->Compute17Next(index + 1);
		if (zh_g.go_back) return;
	}
	for (int idig = 0; idig < 9; idig++) {
		if (idig == digit)continue;
		if (ISFALSEON && zhgxn.nua) continue;
		if (FD[idig][0].On(xcell)) {
			ZHOU* mynext = (this + 1);
			*mynext = *this;
			if (cell >= 54)mynext->ISFALSEON++;
			mynext->SetaCom(idig, cell, xcell);
			mynext->Compute17Next(index + 1);
			if (zh_g.go_back) return;
		}
	}
}



void G17B::Out17(uint32_t bfb3) {
	char ws[82];
	strcpy(ws, empty_puzzle);
	for (int i = 0; i < nclues; i++) {
		int cell = tclues[i];
		ws[cell] = grid0[cell] + '1';
	}
	for (int i = 0, bit = 1; i < 27; i++, bit <<= 1)if (bfb3 & bit)
		ws[i+54]= grid0[i+54] + '1';
	cout << ws << "\t\t one sol   "  << endl;
	fout1 << ws << ";" << genb12.nb12 / 64 << ";" << genb12.i1t16 << ";" << genb12.i2t16 << endl;
	a_17_found_here++;

}
