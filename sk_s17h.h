#define stack1_54 07007007007007007
struct OPCOMMAND {// decoding command line option for this rpocess
	// processing options 
	int opcode;
	int tbn, bx3;// 
	int t18, p1, p2, p2b,//17 of 18 clues, pass or 2 (2a or 2b)
		out_one,// limit output to one per band 3 .bfx[2] & 1
		out_entry; //output of the entry file for test DLL .bfx[2] & 2
	// bfx[2] & 8 special use b2_is as limit b3
	int b1;//band 1 in process 
	int b2,b2_is ;//bands b2  forced
	char* b2start;
	int first, last;
	int ton;//test on and test level
	uint64_t f0, f3, f4, f7,f10; // filters p_cpt2g [3] [4) [7]
	int upto3, upto4; // active below f3 below f4
	void SetUp(int opcod,int k = 0,int p=1) {// init known or not
		memset(this, 0, sizeof * this);
		opcode = opcod;
		tbn = sgo.vx[10];
		bx3 = sgo.vx[11];
		f0= sgo.vx[12];
		if (sgo.bfx[0] & 1)t18 = 1;
		if (sgo.bfx[0] & 6) {// pass 1 2a 2b
			p2 = 1;
			if (sgo.bfx[0] & 4) p2b = 1;
			if (p2b && t18) { p2b = 0; }
		}
		else p1 = 1;
		b1 = sgo.vx[0];
		first = sgo.vx[2];
		last = sgo.vx[3];
		if (last == 0) if (last == 0) { last = first; sgo.vx[3] = first; }  // default if no -v3- given
		b2_is = sgo.vx[4];
		b2 = sgo.vx[5];
		if (sgo.s_strings[0])	if(strlen(sgo.s_strings[0]))
			b2start = sgo.s_strings[0];
		ton= sgo.vx[1];

		f3 = sgo.vx[6];	f4 = sgo.vx[7];		
		f7 = sgo.vx[8]; f10 = sgo.vx[9];

		if (sgo.bfx[1] & 1)upto3 = 1;		
		if (sgo.bfx[1] & 2)upto4 = 1;

		if (sgo.bfx[2] & 1) out_one = 1;
		if (sgo.bfx[2] & 2) out_entry = 1;
		// sgo.bfx[3] is for partial process 

		if (p) {
			cout << Char9out(sgo.bfx[0]) << " sgo.bfx[0 " << endl;
			cout << "standard processing commands_______________" << endl;
			if(t18) cout <<"\t\tsearch 18 clues via -b0-x."<<endl;
			else cout << "\t\tsearch 17 clues via -b0-x." << endl;
			if(p1)cout << "\t\tpass1 via -b0-.x." << endl;
			if (p2)cout << "\t\tpass2 via -b0-.x." << endl;
			if (p2b)cout << "\t\tpass2b via -b0-..x." << endl;
			cout << sgo.vx[0] << " b1  -v0- band 0_415" << endl;
			cout << sgo.vx[2] << " first  -v2- first slice to process" << endl;
			cout << sgo.vx[3] << " last  -v3- last slice to process must be >= vx[2]" << endl;
			if (out_one) cout << " max one out per band 3 sgo.bfx[2] & 1 " << endl;
			if (out_entry)  cout << " file1 contains attached solution grids" << endl;
			cout << "debugging commands___________________" << endl;
			if(b2)		cout << b2 << " b2 -v5- filter band 2 index" << endl;
			if (bx3<416)		cout << bx3 << " b3 index" << endl;
			if (b2start)	cout << b2start << " filter band 2 start" << endl;

			if (ton)cout << ton << "  test on  -v1- verbose mode " << endl;
			if (f0)cout << f0 << "  f0  -vz- diag skip pairs b1b2" << endl;
			if (f3)cout << f3 << "  f3  -v6- diag filter 3 clues [3]" << endl;
			if (f4)cout << f4 << "  f4  -v7- diag filter 6 clues [6]" << endl;
			if (f7)cout << f7 << "  f7  -v8- diag filter go full [7]" << endl;
			if (upto3)cout << "upto debugging [3]  sgo.bfx[1] & 1 " << endl;
			if (upto4)cout << "upto debugging [4]  sgo.bfx[1] & 2 " << endl;

		}
	}
}op;

#define UACNBLOCS 15
//___ expand uas in bands 1+2
struct SPB03A {// spots 6 first clues
	BF128 v;	
	uint64_t  possible_cells, all_previous_cells, active_cells;
};
struct SPB03B {// spots 6 first clues
	BF128 v[UACNBLOCS];
	uint64_t  possible_cells, all_previous_cells, active_cells;
};

// expand uas in band 3
struct SP3 {
	uint32_t  possible_cells, all_previous, active,
		indtw3;
};

//___ permanent data for guas2 guas3
struct SGUA2 {// 81 possible UA2 sockets
	// permanent data
	//uint64_t* tua;
	int col1, col2;// columns of the socket
	int i_81; // index 0_80 for this 
	int i9;// complementary column in minirow
	int id1, id2; // index of digits in gang 27 
	// Current band1+2 data
	int digs, dig1, dig2;// depending on gang27 status
	int valid, // valid if guas 
		validuas,// gua2s found
		used;// if needed in bands3
	int gangcols[9];// revised gangster

}tsgua2[81];
struct SGUA3 {// 81 possible UA3 sockets
	// permanent data
	//uint64_t* tua;// , killer;
	int col1;// first columns 0-9 
	int i_81, stack;// , iguan; // index 0_80 for this 
	int id1, id2, id3; // index of digits in gang 27 
	// Current band1+2 data
	int  dig1, dig2, dig3, digs;// depending on gang27 status
	int valid, // valid if guas 
		validuas,// gua2s found
		used;// if needed in bands3
}tsgua3[81];


#define UA12SIZE 3840
#define UA12BLOCS 30

//_____ uas bands 1+2
struct TUASB12 {//  initial set of uas bands 1+2 
	uint64_t tua[UA12SIZE]; // 30x128
	uint32_t nua,  tdigs[UA12SIZE],ndigs[UA12SIZE];
	inline void AddInit(
		int64_t ua, uint32_t digs,  uint32_t nd) {
		tdigs[nua] = digs;
		ndigs[nua] = nd;
		tua[nua++] = ua;
	}
}tuasb12;
struct T54B12 {//   uas bands 1+2 in 54 mode
	struct TUVECT {//  128 uas and vectors
		BF128 v0, vc[54];
		uint64_t t[128];
		void Init() {
			v0.SetAll_0();
			memset(vc, 255, sizeof vc);
		}
		inline void DoAnd(BF128 v, uint64_t & wa) {
			int ir;
			while ((ir = v.getFirst128()) >= 0) {
				v.clearBit(ir);
				wa &= t[ir];
			}
		}
	};
	// working area for "build"
	uint64_t tw[UA12SIZE]; // to check redundancy
	BF128 vsize[25][UA12BLOCS];
	BF128 tvw[UA12BLOCS];

	// initial status after harvest plus fresh uas B12
	TUVECT ta128[UA12BLOCS];// max start 20*128=2560 uas 
	uint32_t na128, nablocs, nta128[UA12BLOCS];
	void InitA() {
		memset(nta128, 0, sizeof nta128);
		na128 = nablocs = 0;
		for (int i = 0; i < UA12BLOCS; i++)ta128[i].Init();
	}
	inline void AddA(uint64_t u) {
		if (na128 >= UA12BLOCS * 128) return;
		register uint32_t bloc = na128 >> 7, ir = na128 - (bloc << 7);
		na128++; nablocs = bloc; nta128[bloc]++;
		ta128[bloc].v0.setBit(ir);
		BF128* myvc = ta128[bloc].vc;
		register uint64_t R = u;
		ta128[bloc].t[ir] = R;
		R &= BIT_SET_54;// clear extra bits
		uint32_t cell;
		while (bitscanforward64(cell, R)) {
			R ^= (uint64_t)1 << cell; //clear bit
			myvc[cell].clearBit(ir);
		}
	}
	void Build_ta128(uint64_t* t, uint32_t n);

#define UABNBLOCS 20
	// status after 3 clues (B)
	TUVECT tb128[UABNBLOCS];// max start 10*128=1280 uas 
	uint32_t nb128, nbblocs, ntb128[UABNBLOCS];
	void InitB() {
		memset(ntb128, 0, sizeof ntb128);
		nb128 = nbblocs = 0;
		for (int i = 0; i < UABNBLOCS; i++)tb128[i].Init();
	}
	inline void AddB(uint64_t u) {
		if (nb128 >= UABNBLOCS*128) return;
		register uint32_t bloc = nb128 >> 7,
			ir = nb128 - 128 * bloc;
		nb128++; nbblocs = bloc; ntb128[bloc]++;
		tb128[bloc].v0.setBit(ir);
		BF128* myvc = tb128[bloc].vc;
		register uint64_t R = u;
		tb128[bloc].t[ir] = R;
		register uint32_t cell;
		while (bitscanforward64(cell, R)) {
			R ^= (uint64_t)1 << cell; //clear bit
			myvc[cell].clearBit(ir);
		}
	}
	int Build_tb128( SPB03A& s);

//   #define UACNBLOCS 15 see befor SPB 
	// status after 6 clues (C)
	TUVECT tc128[UACNBLOCS];// max start 10*128=1280 uas 
	uint64_t tandc;
	uint32_t nc128, ncblocs, ntc128[UACNBLOCS];
	void InitC() {
		memset(ntc128, 0, sizeof ntc128);
		nc128 = ncblocs = 0;
		for (int i = 0; i < UACNBLOCS; i++)tc128[i].Init();
		tandc = ~0;
	}
	inline void AddC(uint64_t u) {
		if (nc128 >= UACNBLOCS * 128) return;
		register uint32_t bloc = nc128 >> 7,
			ir = nc128 - 128 * bloc;
		nc128++; ncblocs = bloc; ntc128[bloc]++;
		tc128[bloc].v0.setBit(ir);
		BF128* myvc = tc128[bloc].vc;
		register uint64_t R = u;
		tc128[bloc].t[ir] = R;
		tandc &= R;
		register uint32_t cell;
		while (bitscanforward64(cell, R)) {
			R ^= (uint64_t)1 << cell; //clear bit
			myvc[cell].clearBit(ir);
		}
	}
	int Build_tc128( SPB03A& s6);// after 6 clues
	inline int IsNotRedundant(uint64_t u) {
		register uint64_t nu = ~u;
		for (uint32_t i = 0; i < ntc128[0]; i++)
			if (!(tc128[0].t[i] & nu)) return 0;
		return 1;
	}

#define UADNBLOCS 2
	// nearly no chance to get a 12 with more than 64 uas
	// status after 9 clues (D)
	TUVECT td128[UADNBLOCS];// max start 10*128=1280 uas 
	uint64_t tandd;
	uint32_t nd128, ndblocs, ntd128[UADNBLOCS];
	void InitD() {
		memset(ntd128, 0, sizeof ntd128);
		nd128 = ndblocs = 0;
		for (int i = 0; i < UADNBLOCS; i++)td128[i].Init();
		tandd = ~0;
	}
	inline void AddD(uint64_t u) {
		if (nd128 >= UADNBLOCS * 128) return;
		register uint32_t bloc = nd128 >> 7,
			ir = nd128 - ( bloc<<7);
		nd128++; ndblocs = bloc; ntd128[bloc]++;
		td128[bloc].v0.setBit(ir);
		BF128* myvd = td128[bloc].vc;
		register uint64_t R = u;
		td128[bloc].t[ir] = R;
		tandd &= R;
		register uint32_t cell;
		while (bitscanforward64(cell, R)) {
			R ^= (uint64_t)1 << cell; //clear bit
			myvd[cell].clearBit(ir);
		}
	}
	int Build_td128(SPB03A& s9);
	inline int IsNotRedundantD(uint64_t u) {
		register uint64_t nu = ~u;
		for (uint32_t i = 0; i < ntd128[0]; i++)
			if (!(td128[0].t[i] & nu)) return 0;
		return 1;
	}

}t54b12;

//____ guas (uas with band3)
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
			if (nua < 128)tua[nua++] = ua;
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
		void SortClean(uint32_t lim=40) {
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
					if (nua > lim) break;
				}
				if (nua > lim) break;
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
	void SortClean(uint32_t lim = 40) {
		for (int i = 0; i < 81; i++)
			if (tg2[i].nua > 1) tg2[i].SortClean(lim);
	}

	int CutG2x(int lim) {
		int n = 0;
		for (int i = 0; i < 81; i++) {
			if ((int)tg2[i].nua > lim)tg2[i].nua = lim;
			n += tg2[i].nua;
		}
		return n;
	}
	int CutG3x(int lim) {
		int n = 0;
		for (int i = 0; i < 81; i++) {
			if ((int)tg3[i].nua > lim)tg3[i].nua = lim;
			n += tg3[i].nua;
		}
		return n;
	}

}guah;
struct GUAH54N {
	struct Z128 {// for one gua2/gua3
		BF128 	vc[54], vi[27],// vectore/cell /i81
			v0, //active items in vector
			v6, v9; // v0 after clues
		uint64_t killer;
		uint32_t mode2_3, i81, n;
		void Init(uint32_t m2_3, uint32_t ix) {
			memset(vc, 255, sizeof vc);
			memset(&v0, 0, sizeof v0);
			mode2_3 = m2_3;
			i81 = ix;
			n = 0;
			killer = ~0;
		}
		void Enter(uint64_t u) {
			if (n >= 128) return;
			uint32_t nr = n++;
			killer &= u;
			v0.setBit(nr);
			v6.clearBit(nr);// free if fresh ua
			v9.clearBit(nr);// free if fresh ua
			register uint32_t cell;
			register uint64_t U = u;
			while (bitscanforward64(cell, U)) {
				U ^= (uint64_t)1 << cell;
				vc[cell].clearBit(nr);
			}
		}

		int Redundant(uint64_t u) {
			for (uint32_t i = 0; i < 128; i++) {
				if (v0.Off(i))return 0;
				register uint64_t uc = 0, bit = 1;
				for (int j = 0; j < 54; j++,bit<<=1)
					if (vc[j].Off(i)) uc|=bit;
				if(u==uc )return 1;
			}
			return 0;
		}

	}zz[162];
	struct BLOC0 {// only one ua maxi 16
		uint64_t ua, id;
	}bloc0[40];
	uint32_t nbl0, nzzused, nzzg2;
	uint32_t tind6[162], ntind6;
	uint32_t tind9[162], ntind9;
	uint32_t tg2[81], ntg2,tmg2[20], ntmg2;
	uint32_t tg3[81], ntg3,tmg3[20], ntmg3;
	BF128 g2, g3,g2_6,g3_6,g2_9,g3_9;
	int indg2[81], indg3[81];

	//____________________________________
	void Init() {	
		nbl0 = nzzused =0;	
		memset(indg2, 255, sizeof indg2);
		memset(indg3, 255, sizeof indg3);
	}
	inline void InitCom() { ntmg2 = ntmg3 = 0; }
	void Build();
	void Build6( uint32_t* tc) {
		g2_6.SetAll_0(); g3_6 = g2_6;
		ntind6 = 0;
		for (uint32_t i = 0; i < nzzused; i++) {
			Z128& myz = zz[i];
			//BF128 v = myz.v0;
			BF128 v= myz.vc[tc[0]];
			for (int j = 1; j < 6; j++)
				v &= myz.vc[tc[j]];
			myz.v6 = v; myz.v9 = v;
			if ((v & myz.v0).isNotEmpty()) {
				tind6[ntind6++] = i;
				if (myz.mode2_3) g3_6.setBit(myz.i81);
				else g2_6.setBit(myz.i81);
			}
		}
	}
	void Build9(uint32_t* tc) {
		g2_9.SetAll_0(); g3_9 = g2_9;
		ntind9 = 0;
		for (uint32_t it = 0; it < ntind6; it++){
			Z128& myz = zz[tind6[it]];
			BF128 v = myz.v6;
			for (int j = 0; j < 3; j++)
				v &= myz.vc[tc[j]];
			myz.v9 = v;
			if ((v & myz.v0).isNotEmpty()) {
				tind9[ntind9++] = tind6[it];
				if (myz.mode2_3) g3_9.setBit(myz.i81);
				else g2_9.setBit(myz.i81);
			}
		}
	}
	inline void GetG2G3A(uint64_t bf) {
		g2.SetAll_0(); g3 = g2;
		ntg2 = ntg3 = 0;
		// see bloc 0
		for (uint32_t i = 0; i < nbl0; i++) {
			if (!(bf & bloc0[i].ua)) {
				register uint32_t ix = (int)bloc0[i].id;
				if (ix > 100) {// this is a g3
					ix -= 100;
					g3.setBit(ix);
					tg3[ntg3++] = ix;
				}
				else {
					g2.setBit(ix);
					tg2[ntg2++] = ix;
				}
			}
		}
	}
	void GetG2G3_7(uint64_t bf, uint32_t cell1) {
		GetG2G3A(bf);
		{
			for (uint32_t it = 0; it < ntind6; it++) {
				Z128& myz = zz[tind6[it]];
				if (bf & myz.killer) continue;
				BF128 v = (myz.v6 & myz.v0) & myz.vc[cell1];
				if (v.isNotEmpty()) {
					register uint32_t ix = myz.i81;
					if (myz.mode2_3) {// this is a g3
						g3.setBit(ix);
						tg3[ntg3++] = ix;
					}
					else {
						g2.setBit(ix);
						tg2[ntg2++] = ix;
					}
				}
			}
		}
	}
	void GetG2G3_8(uint64_t bf, uint32_t * tc) {
		GetG2G3A(bf);
		{
			for (uint32_t it = 0; it < ntind6; it++) {
				Z128& myz = zz[tind6[it]];
				if (bf & myz.killer) continue;
				BF128 v = (myz.v6 & myz.v0) & 
					(myz.vc[tc[0]]& myz.vc[tc[1]]);
				if (v.isNotEmpty()) {
					register uint32_t ix = myz.i81;
					if (myz.mode2_3) {// this is a g3
						g3.setBit(ix);
						tg3[ntg3++] = ix;
					}
					else {
						g2.setBit(ix);
						tg2[ntg2++] = ix;
					}
				}
			}
		}
	}
	void GetG2G3_9(uint64_t bf) {
		GetG2G3A(bf);
		for (uint32_t it = 0; it < ntind9; it++) {
			Z128& myz = zz[tind9[it]];
			register uint32_t ix = myz.i81;
			if (myz.mode2_3) {// this is a g3
				g3.setBit(ix);
				tg3[ntg3++] = ix;
			}
			else {
				g2.setBit(ix);
				tg2[ntg2++] = ix;
			}
		}
	}
	void GetG2G3_10(uint64_t bf, uint32_t cell1) {
		GetG2G3A(bf);
		{
			for (uint32_t it = 0; it < ntind9; it++) {
				Z128& myz = zz[tind9[it]];
				if (bf & myz.killer) continue;
				BF128 v = (myz.v9 & myz.v0) & myz.vc[cell1];
				if (v.isNotEmpty()) {
					register uint32_t ix = myz.i81;
					if (myz.mode2_3) {// this is a g3
						g3.setBit(ix);
						tg3[ntg3++] = ix;
					}
					else {
						g2.setBit(ix);
						tg2[ntg2++] = ix;
					}
				}
			}
		}
	}
	void GetG2G3_11(uint64_t bf, uint32_t cell1, uint32_t cell2) {
		GetG2G3A(bf);
		{
			for (uint32_t it = 0; it < ntind9; it++) {
				Z128& myz = zz[tind9[it]];
				if (bf & myz.killer) continue;
				BF128 v = (myz.v9 & myz.v0) & (myz.vc[cell1]& myz.vc[cell2]);
				if (v.isNotEmpty()) {
					register uint32_t ix = myz.i81;
					if (myz.mode2_3) {// this is a g3
						g3.setBit(ix);
						tg3[ntg3++] = ix;
					}
					else {
						g2.setBit(ix);
						tg2[ntg2++] = ix;
					}
				}
			}
		}
	}
	void GetG2G3_12(uint64_t bf, uint32_t* tc) {
		GetG2G3A(bf);
		{
			for (uint32_t it = 0; it < ntind9; it++) {
				Z128& myz = zz[tind9[it]];
				if (bf & myz.killer) continue;
				BF128 v = (myz.v9 & myz.v0);
				for (uint32_t j = 0; j < 3; j++)
					v &= myz.vc[tc[j]];
				if (v.isNotEmpty()) {
					register uint32_t ix = myz.i81;
					if (myz.mode2_3) {// this is a g3
						g3.setBit(ix);
						tg3[ntg3++] = ix;
					}
					else {
						g2.setBit(ix);
						tg2[ntg2++] = ix;
					}
				}
			}
		}
	}

	inline void AddA2(uint64_t bf,  int32_t i81,int cc12){
		
		register int iz = indg2[i81];
		if (g2_9.Off(i81)) {	
			g2_9.setBit(i81);
			tind9[ntind9++] = iz;
			tmg2[ntmg2++] = i81;
			if (g2_6.Off(i81)) {
				g2_6.setBit(i81);
				tind6[ntind6++] = iz;
			}
		}
		else {
			uint32_t nx = ntmg2;
			tmg2[ntmg2++] = i81;
			for (uint32_t i = 0; i < nx; i++) {
				if (i81 == tmg2[i]) {
					ntmg2--; // was there
					break;
				}
			}
		}
		int n =  zz[iz].n;
		if (n >= 80 && cc12 > 14) return;
		if (n >= 100 && cc12 > 13) return;
		if(n>120 && n<128)
		cout << "\t\t\tadded g2 " << iz << " " << zz[iz].n << endl;
		zz[iz].Enter(bf);
	}
	inline void AddA3(uint64_t bf, int32_t i81) {
		register int ix = indg3[i81];

		if (g3_9.Off(i81)) {
			g3_9.setBit(i81);
			tind9[ntind9++] = ix;
			tmg3[ntmg3++] = i81;
			if (g3_6.Off(i81)) {
				g3_6.setBit(i81);
				tind6[ntind6++] = ix;
			}
		}
		else {
			uint32_t nx = ntmg3;
			tmg3[ntmg3++] = i81;
			for (uint32_t i = 0; i < nx; i++) {
				if (i81 == tmg3[i]) {
					ntmg3--; // was there
					break;
				}
			}
		}

		int n = zz[ix].n;
		zz[ix].Enter(bf);
	}

	void Add(uint64_t bf, BF128* vc, BF128 & v0,
		char * ix,	int32_t n, int32_t i81) {
		v0.setBit(n);
		ix[n] = i81;
		register uint32_t cell;
		register uint64_t U = bf;
		while (bitscanforward64(cell, U)) {
			U ^= (uint64_t)1 << cell;
			vc[cell].clearBit(n);
		}
	}

}guah54n;

//_____ processing band3
struct XQ {//to build the UAs b3 to expand
	uint32_t t1a; //27 bits field common 2 pairs 
	uint32_t  critbf,fa,fb;
	uint32_t  * t2a, * t2b3,*t2b;
	uint32_t n2a, n2b,n2b3, nb3,bf3p,nmiss,nadded;
	uint32_t nin, nout, nred;
	uint32_t iuas4;
	uint32_t tbuf[60];
	uint32_t tin[400], tout[400],tred[300];
	void Init() { 
		t1a = critbf = 0;
		n2a = n2b = n2b3= bf3p = 0;
		nin=nout=nadded= 0; 
		t2a = tbuf; t2b3 = tbuf + 12; t2b = t2b3 + 24;
	}
	inline void SetFilters () {
		if (!nmiss) {
			if ((fa = t1a)) {// initial safety
				critbf = 0;
				for (uint32_t i = 0; i < n2b3; i++)
					critbf |= t2b3[i];
				for (uint32_t i = 0; i < n2b; i++)
					critbf |= t2b[i];
			}
			fb = critbf;
		}
		else {	fa = 0; fb = BIT_SET_27;}
	}
	int Miss1ToMiss0();
	int MissxToMiss0(uint32_t ubf);
	int Miss0CheckTin();
	int NoRoomToAssign() {
		return(_popcnt32(t1a) >= nb3);
	}
	int NToAssign() {
		return( nb3- _popcnt32(t1a));
	}
	inline void SetFreshCrit() {
		fb = 0;
		for (uint32_t i = 0; i < n2b; i++) {
			fb |= t2b[i];
		}
		critbf = fb;
	}
	uint32_t  AssignMiss0(uint32_t bf);

	uint32_t  AddAssigned(uint32_t bf) {
		nadded++;
		t1a |= bf;
		fa = t1a;
		register uint32_t F = t1a, C = 0,n=n2b;
		n2b = 0;
		for (uint32_t i = 0; i < n; i++) {
			register uint32_t U = t2b[i];
			if (!(U & F)) {
				t2b[n2b++]=U;
				C |= U;
			}
		}
		fb=critbf = C;
		return C;
	}
	inline void Addin(uint32_t bf) { tin[nin++] = bf; }
	inline void Addout(uint32_t bf) { tout[nout++] = bf; }
	inline uint32_t GetAndout() {
		register uint32_t A = tout[0];
		for (uint32_t i = 1; i < nout; i++)
			A&=tout[i];
		return A;
	}
	int AddRedundant(uint32_t bf) {
		register uint32_t nbf = ~bf;
		for (uint32_t i = 0; i < nred; i++) {
			register uint32_t u = tred[i];
			if (!(u & nbf)) return 0;// subset or =
			if (!(bf & ~u)) {
				tred[i]=bf; // new is subset
				return 1;// keep it as valid
			}
		}
		tred[nred++] = bf;
		return 1;
	}
	void BuildMiss0Redundant() {
		nred = n2b3+n2b;
		memcpy(tred, t2b3, n2b3 * sizeof tred[0]);
		memcpy(&tred[n2b3], t2b, n2b * sizeof tred[0]);
		int n = 0;
		for (uint32_t i = 0; i < nin; i++) {
			if (AddRedundant(tin[i]))
				tin[n++] = tin[i];;
		}
		nin = n;
	}
	void BuildMissxRedundant() {
		nred =n2a+ n2b3 + n2b;
		memcpy(tred, t2a, n2a * sizeof tred[0]);
		memcpy(&tred[n2a], t2b3, n2b3 * sizeof tred[0]);
		memcpy(&tred[n2a+n2b3], t2b, n2b * sizeof tred[0]);
		int n = 0;
		for (uint32_t i = 0; i < nin; i++) {
			if (AddRedundant(tin[i]))
				tin[n++] = tin[i];;
		}
		nin = n;
	}
	inline int SubCritMore(uint32_t F) {
		if (!F) return 0;
		register uint32_t ff = F, ret = 0;
		for (uint32_t i = 0; i < n2a; i += 2) {
			register uint32_t a = t2a[i], b = t2a[i + 1],
				mask = a | b, f = mask & F;
			if (f) {
				ff &= ~(mask & F);
				int cc = _popcnt32(f);
				if (cc > 1)ret += _popcnt32(f) - 1;
				else if (f &= ~(a & b)) ret++;
				//f &= ~(a & b); // clues extra mincount
				//ret += _popcnt32(f)-1;
			}
		}
		if (ret > 1 || !ff) return ret; 
		for (uint32_t i = 0; i < n2b3; i += 3) {
			register uint32_t mask = t2b3[i] | t2b3[i + 1]| t2b3[i + 2],
				f = mask & F;
			if (f) {
				if (f == mask) ret++;
				ff &= ~(mask & F);	
			}
		}
		if (ret > 1||!ff) return ret;// already to many
		// others are disjoints
		for (uint32_t i = 0; i < n2b; i++) {
			register uint32_t a = t2b[i] & ff;
			if (a) 	ret += (_popcnt32(a) - 1);
		}
		return ret;
	}
	void SubCritMoreDo(uint32_t F) {// 1 more to process
		// clean all hits
		register uint32_t nn = n2a;
		n2a = 0;
		for (uint32_t i = 0; i < nn; i += 2) {
			register uint32_t a = t2a[i], b = t2a[i + 1];
			if (!((a | b ) & F)) {
				t2a[n2a++] = a; t2a[n2a++] = b;
			}
			else {
				if (!(a & F))t2b[n2b++] = a;
				if (!(b & F))t2b[n2b++] = b;
				t1a &= ~(a | b);
			}
		}
		nn = n2b3; n2b3 = 0; // Is it 3 pairs one minirow
		for (uint32_t i = 0; i < nn; i += 3) {
		 	register uint32_t a= t2b3[i],
				b= t2b3[i + 1],c= t2b3[i + 2];
			if (!((a | b | c) & F)) {
				t2b3[n2b3++] = a;
				t2b3[n2b3++] = b;
				t2b3[n2b3++] =c;
			}
			else {
				if (!(a & F))t2b[n2b++] = a;
				if (!(b & F))t2b[n2b++] = b;
				if (!(c & F))t2b[n2b++] = c;
			}	
		}
		// others are disjoints more if 2 clues
		nn = n2b;	n2b = 0;
		for (uint32_t i = 0; i < nn; i++) {
			register uint32_t a = t2b[i] ;
			if (!(a & F))t2b[n2b++] = a;
		}

	}
	inline uint32_t SubCritActive(uint32_t F) {
		register uint32_t ret = 0;
		for (uint32_t i = 0; i < n2a; i += 2) {
			register uint32_t a = t2a[i], b = t2a[i + 1],
				mask = a | b, f = mask & F;
			if (!f) ret |= mask;
		}
		for (uint32_t i = 0; i < n2b3; i ++)  
			if(!(F& t2b3[i]))ret |= t2b3[i];
		for (uint32_t i = 0; i < n2b; i++)
			if (!(F & t2b[i]))ret |= t2b[i];		
		return ret;
	}
	inline uint32_t NonHitOut(uint32_t F) {
		register uint32_t ret = 0;
		for (uint32_t i = 0; i < nout; i ++) 
			if (!(F&tout[i]) )return 1;
		return 0;
	}
	int IsNotRedundantOut(uint32_t bf) {
		register uint32_t nbf = ~bf;
		for (uint32_t i = 0; i < nout; i++) {
			register uint32_t u = tout[i];
			if (!(u & nbf)) return 0;// subset or =			
		}
		return 1;
	}
	void BuildMiss0Out() {
		nout = n2b3;
		nred = 0;
		memcpy(tout, t2b3, n2b3 * sizeof tred[0]);
		// pickup all 2 store more 
		for (uint32_t i = 0; i < n2b; i++) {
			register uint32_t u = t2b[i];
			if (_popcnt32(u) > 2) {
				tred[nred++] = u;
				continue;
			}
			tout[nout++] = u;
		}
		for (uint32_t i = 0; i < nin; i++) {
			register uint32_t u = tin[i];
			if (_popcnt32(u) > 2) {
				tred[nred++] = u;
				continue;
			}
			if(IsNotRedundantOut(u))
				tout[nout++] = u;
		}
		for (uint32_t i = 0; i < nred; i++) {
			register uint32_t u = tred[i];
			if (IsNotRedundantOut(u))
				tout[nout++] = u;
		}

	}
	void BuildMissxOut() {
		register uint32_t* tt = tout + nout;
		memcpy(tt, t2a, n2a * sizeof tred[0]);
		tt += n2a;
		memcpy(tt, t2b3, n2b3 * sizeof tred[0]);
		tt += n2b3;
		memcpy(tt, t2b, n2b * sizeof tred[0]);
		nout =(int) (tt - tout) + n2b;
		for (uint32_t i = 0; i < nin; i++) {
			register uint32_t u = tin[i];
			if (IsNotRedundantOut(u))
				tout[nout++] = u;
		}
		nin = 0;
	}


	void BuildCheckRedundant() {
		nred = n2b;
		if (nmiss) {
			nred += n2a;
			memcpy(&t2b[n2b], t2a, n2a * sizeof t2b[0]);
		}
	}
	uint64_t Isoutsize2();
	uint64_t Isoutsize3();
	uint64_t Isoutsize4();
	int NoDisjoint() {
		register uint32_t r1 = tout[0] , r2 ;
		for (uint32_t i = 1; i < nout; i++) {
			register uint32_t u = tout[i];
			if (r2 |= r1 & u)return 1; 
				r1 |= u;
		}
		return 0;
	}
	int Min1_4Disjoint();
	void CleanIn();
	void CleanOut();
	void CleanOut(uint32_t F, uint32_t A);
}xq;

// standard  bands  
struct STD_B416 {
	char band[28];// band in char mode
	int i416,// id 0-415 in minlex mode of the band
		map[27],// mapping cells from minlex to solution grid
		band0[27], // band in 0-8 integer mode
		gangster[9], // bit field digits per column
		dpat[9],// solution pat per digit
		dpmi[9],// initial map gangster per digit
		dband;
	uint32_t tua[82], nua;//  maximum 81 morphed uas 
	uint32_t fd_sols[2][9];//start puzzle/ solution
	void Initstd();// first band initial (once)
	void GetBandTable(int i416e);//pick up band 1 from table
	void InitG12(int i416e);// end init band 1 
	void InitBand2_3(int i16, char* ze, BANDMINLEX::PERM& p
		, int iband = 1);

	void SetGangster();
	inline void GetUAs() {
		nua = t16_nua[i416];
		memcpy(tua, &t16_UAs[t16_indua[i416]], 4 * nua);
	}
	void MorphUas();
	void InitC10(int i);// known mode

};
struct STD_B3 :STD_B416 {// data specific to bands 3
	// permanent gangster information
	struct G {
		BF128 gsocket2, gsocket3;// active i81 mode 81 bits
		BF128 gsocket4, gsocket6;//g2 with 4 cells in band 3
		BF128 gsocket2all;
		int pat2[81], pat3[81]; // storing ua bitfields
		int pat2_27[27]; // storing ua bitfields
		int ua2_imini[81], ua3_imini[81], ua2bit27[81];
	}g;
	char isg2[81],isg3[81];
	uint32_t  mg2, mg3,minc,ming2_1, ming2_2, ming2_3;
	int minirows_bf[9];
	int triplet_perms[9][2][3];
	uint32_t i_27_to_81[27], i_9_to_81[9]; //band i81 for guas guas3
	uint32_t i_81_to_27[81]; //band i81 for guas guas3
	uint32_t  poutdone, aigskip,
		tg2_4[50], ntg2_4, tg2_6[50], ntg2_6;
	//_______________________
	void InitBand3(int i16, char* ze, BANDMINLEX::PERM& p);
	void GoA();// start band 3  
	void GoAg23();// fresh g2 g3 to consider 
	void GoB0();// band3 miss0 at start
	void GoC0(uint32_t bf, uint32_t a);// band3 miss0 before guam guamm
	void GoC0F(uint32_t bf);//same nb3 filled
	void GoB1();//  band 3 miss1 at start
	void GoB1toMiss0(uint32_t bf);//  band 3 miss1 at start
	void GoBMore1();//  band 3 miss>1 at start
	void GoBMoretoMiss0(uint64_t ubf);
	void GoBMoretoMiss1(uint64_t ubf);


	uint32_t Get2d(int d1, int d2) {
		return fd_sols[0][d1] | fd_sols[0][d2];
	}
	inline int GetI81_2(int bf) {
		for (uint32_t i = 0; i < 27; i++) {
			register uint32_t i81 = i_27_to_81[i];
			if (g.pat2[i81] == bf) return i81;
		}
		return -1;
	}
	inline int GetI81_3(int bf) {
		for (uint32_t i = 0; i < 9; i++) {
			register uint32_t i81 = i_9_to_81[i];
			if (g.pat3[i81] == bf) return i81;
		}
		return -1;
	}
	inline int GetI81_x(int bf) {
		for (uint32_t i = 0; i < 81; i++) {
			if (g.pat2[i] == bf) return i;
		}
		return -1;
	}


	void Pat_to_digitsx(int bf, int* tcells, int* tdigits, int& nt) {
		register int cell;
		nt = 0;
		while (bitscanforward(cell, bf)) {
			bf ^= 1 << cell;
			tcells[nt] = cell;
			tdigits[nt++] = band0[cell];
		}
	}
	int Is_Pat_For_Mex(int* tcells, int* tdigits, int nt) {
		for (int i = 0; i < nt; i++) {
			if (band0[tcells[nt]] != tdigits[nt]) return 0;
		}
		return 1;
	}


#define GM_NB4 50
#define GM_NBM 50

	struct GUM64 {
		uint64_t v0,  vc[55];//vc[54] is for 6 clues
		uint32_t tb3[64], n;
		void Init() {
			n = 0;		v0 = 0;
			memset(vc, 255, sizeof vc);
			vc[54] = 0;
		}
		int Add(uint64_t bf, uint32_t patb3) {
			if (n >= 64) return 0;//safety
			tb3[n] = patb3;
			uint64_t bit = (uint64_t)1 << n++;
			v0 |= bit; vc[54] |= bit;
			register uint64_t B = bf, nbit = ~bit;
			int cell;
			while (bitscanforward64(cell, B)) {
				B ^= (uint64_t)1 << cell;
				vc[cell] &= nbit;
			}
			return 1;
		}
		inline void Set6(uint32_t* t) {
			register uint64_t v = v0;
			for (int i = 0; i < 6; i++) v &= vc[t[i]];
			vc[54] = v;
		}
		inline uint64_t Getv(uint32_t* t, int nt) {
			register uint64_t v = vc[54];
			for (int i = 0; i < nt; i++) v &= vc[t[i]];
			return v;
		}

		void Dump(int i0) {
			for (uint32_t i = 0; i < n; i++) {
				uint64_t bit = (uint64_t)1 << i;
				cout <<Char27out(tb3[i]) << "\t";
				for (int j = 0; j < 54; j++)
					if (vc[j] & bit) cout << ".";
					else cout << "1";
				cout << "  " << i + i0 << endl;
			}
		}

		uint64_t bf12;// mode 54
		uint32_t bf3;
		inline int Count() {
			return ((int)_popcnt64(bf12) + _popcnt32(bf3));
		}
		void Print() {
			cout << Char54out(bf12) << "\t";
			cout << Char27out(bf3)
				<< " " << Count() << endl;
		}
	}tgm64[GM_NB4], tgm64m[GM_NBM];
	uint32_t nbgm,  nbbgm, nbgmm, nbbgmm;
	void InitTg() {
		nbgm = nbbgm =  0;
		nbgmm = nbbgmm = 0;
		for (int i = 0; i < GM_NB4; i++) tgm64[i].Init();
		for (int i = 0; i < GM_NBM; i++) tgm64m[i].Init();
	}
	inline void Addm4(BF128 &w) {// entry mode 3x
		if (nbgm >= 64 * GM_NB4) return;
		nbbgm = nbgm++ >> 6; // new last bloc 
		register uint64_t U = w.bf.u64[0];
		U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54
		tgm64[nbbgm].Add(U, w.bf.u32[2]);
	}
	inline void Addmm(BF128& w) {// entry mode 3x
		if (nbgmm >= 64 * GM_NBM) return;
		nbbgmm = nbgmm++ >> 6; // new last bloc 
		register uint64_t U = w.bf.u64[0];
		U = (U & BIT_SET_27) | ((U & BIT_SET_B2) >> 5);// mode 54
		tgm64m[nbbgmm].Add(U, w.bf.u32[2]);
	}

	inline void Set6clues(uint32_t * tc) {
		for (uint32_t i = 0; i <= nbbgm; i++) 
			tgm64[i].Set6(tc);
		for (uint32_t i = 0; i <= nbbgmm; i++) 
				tgm64m[i].Set6(tc);
	}

};


//____ entry builder

struct GEN_BANDES_12 {// encapsulating global data 
	STD_B3 bands3[512];
	int nband3;
	int sgchecked[512][81], nsgchecked;// used in pass1 band3 <=
	int modeb12, go_back, diagmore, diagbug, ip20,
		it16, it16_2, it16_3, imin16_1, imin16_2, imin16_3;
	int i1t16, i2t16, i3t16	,// index 416 ordered in increasing size of valid clues 3
		ixmin1, ixmin2, ixmin3, maxnb3;
	char zsol[82], rband2[28];
	int grid0[81], gw[81] , tc[6], ntc;
	int gcheck[82], ib1check, ib2check, ib3check,ibasecheck;
	uint64_t   nb12;
	BANDMINLEX::PERM t_auto_b1[108], // maxi is 107excluding start position
		t_auto_b1b2[108], t_auto_b2b1[108],
		pband2, pband3,  pcheck2 , pcheck3;
	int n_auto_b1, n_auto_b1b2, n_auto_b2b1;

	int cold[9], coldf[9], rowd[6], boxd[6], rowdb3[3], boxdb3[3]; //free digits 
	//_________________ gangster 
	int gangcols[9];// 9 bits fields band3 for each colum (see valid band2)
	int gangb12[9]; // digit bf bands 12 per column

	//================================= functions
	void GetStartB2(int i); // one of the 20 starts 
	void Start(int mode = 0);
	void NewBand1(int iw);
	int F17Novalid1_2();
	int Band2Check();
	int Band3Check();
	int Band2_3CheckNoauto();
	void Find_band2B();
	int ValidBand2();
	void ValidInitGang();
	void Find_band3B(int m10 = 1);
	void Find_band3B_pass1B(int m10 = 1);
	void OutEntry();
	void F3B_See();
	inline void F3B_See_18();// one NED return 1 if equal not loaded
	void F3B_See_Com();// one NED b3 first   
	int F3B_See_Com_FilterDiag(); 
	void  F3B_See_Com_GetCFX();

	//============= loops control for UAs 5;6;7 digits collection (collect more=
	int iband, ibox, iminirow, ibox2, iminirow2, pat1, pat2, ncells;
	int tcells[6], tcols[6];
	int bcols[2][9], mycols[9], myfloors;
	uint64_t mybf;
};


// main process 

struct G17B {// hosting the search in 6 6 5 mode combining bands solutions
	G17B();// initial tasks all commands
	int b3lim,	 aigstop, aigstopxy,nb3_not_found,
		npuz, a_17_found_here ;
	int ng2,ng3;
	int grid0[81];

	//____gangsters, brute force,sockets setup
	//_________________ gangster 
	int gang[9][3]; // gangcols expanded (buildgang ) 3 digits
	int* gang27; // redefines gang[9][3] as 27 integer
	int   gang_digits_cols[9][3];// active cols for a given digit

	//______sockets common to  all bands 3  
	BF128 gsock2, gsock3;
	
	//============================ b12 no more uas to test
	BF128 dvect[2];
	uint64_t ua_ret7p, myb12, myandall,	myac,
		anduab12,build9done, tcluesxpdone,
		clean_valid_done,critical_done;

	uint64_t bf_cl3, bf_cl6, bf_cl9, bf_clc;
	uint64_t ac_cl3, ac_cl6, ac_cl9, ac_clc;
	uint32_t tcluesxp[20],// to access guam
		tclx[10],// to expan 7_+
		tc_1_3[3], tc_1_6[6], tc_4_6[3], tc_7_9[5],// 5 for expand 7 11 
		tc_10_12[3];
	int ncluesxp, ntc_6_9;
	inline void Set3(SPB03A& s3) {
		bf_cl3 = s3.all_previous_cells;
		ac_cl3 = s3.active_cells;
		memcpy(tc_1_6, tc_1_3,  sizeof tc_1_3);
	}
	inline void Set6(SPB03A& s6) {
		bf_cl6 = s6.all_previous_cells;
		ac_cl6 = s6.active_cells;
		memcpy(&tc_1_6[3], tc_4_6, sizeof tc_4_6);
	}
	inline void Set9(SPB03A& s9) {
		bf_cl9 = s9.all_previous_cells;
		ac_cl9 = s9.active_cells;
		//memcpy(&tc_1_6[3], tc_4_6, sizeof tc_4_6);
	}
	//============  go band3
	int nclgo, nmiss;
	int  ncluesb12, ncluesb3, mincluesb3;
	uint32_t anduab3;// b3 expand
	STD_B3* myband3;

	uint32_t  t3infield, t3outseen;
	uint32_t t3b[200], nt3b, t3c[100], nt3c;// after 3/6 clues
	uint32_t t3more[200], nt3more;
	uint32_t ntoassb3;

	void Dumpt3b() {
		cout << "t3b n=" << nt3b << endl;
		for (uint32_t i = 0; i < nt3b; i++)
			cout << Char27out(t3b[i]) << endl;
	}

	uint32_t t3[1000], nt3,		t3_2[1000], nt3_2,
		uasb3_1[1000], uasb3_2[1000], uas_in[1000],
		nuasb3_1, nuasb3_2, nuas_in, b3_andout;

	void BuildSortT3b();
	/*
	inline uint32_t GetAndT3o() {
		uint32_t wa = BIT_SET_27;
		for (uint32_t i = 0; i < nt3o; i++)
			wa &= t3o[i];
		return wa;
	}	*/

	
	//=====================process for a new band 2 / set of bands 3
	void Start();// standard entry
	void StartInit();// initial task gangster set up 
	void UaCollector();
	inline void Adduab12(uint32_t digs, uint32_t nd);
	void FirstUasCollect();
	void SecondUasCollect();
	void UasCollect4box();
	void UasCollect6_7();

	void StartAfterUasHarvest();
	//inline int BuildGua(BF128& w);
	//inline void BuildGua(BF128& w, int cc);
	void Guas2Collect();
	void Guas2CollectG2();
	void Guas2CollectG3();
	void Guas2CollectG3_4d();
	void Guas2CollectG3_5d();

	int IsValid7pbf(uint64_t bf);

	void Expand_03();
	void Expand_46(SPB03A& s3);
	void Expand_7_9(SPB03A& s6);
	void EndExpand_7_9();
	void Expand_10_11_17(SPB03A& s9);
	void Expand_10_11_18(SPB03A& s9);
	void Valid9(SPB03A& s9);

	void Expand_10_12(SPB03A& s9);
	void EndExpand_10_12();

	inline void DoBuild9() {
		if (!build9done) {
			build9done = 1;
			guah54n.Build9(tc_7_9);
		}
	}
	void Do7x() {
		if (!tcluesxpdone) {
			tcluesxpdone = 1;
			register uint64_t F = myb12 & ~bf_cl6; // cells 6 to x
			ncluesxp = 0;
			register int cell;// build table of cells 
			while (bitscanforward64(cell, F)) {
				F ^= (uint64_t)1 << cell;
				tcluesxp[ncluesxp++] = cell;
			}
		}
	}

	void GoCallB3Com();

	uint32_t IsValidB3(uint32_t bf);

	// option no table of clues, no live count per band stack


	int IsValid_myb12();

	void Go_8_10();
	void Go_7_10();
	void Go_x_10();

	void Go_8_11_18();
	void Go_7_11_18();
	void Go_x_11_18();


	void Go_8_12();
	void Go_7_12();
	void Go_x_12();

	void Go_8_11_17();
	void Go_7_11_17();
	void Go_x_11_17();

	inline int VerifyValidb3() {
		if (clean_valid_done == 2) return 1;
		if (clean_valid_done) return 0;
		clean_valid_done = 1;
		if (!IsValid_myb12()) return 0;
		clean_valid_done = 2;
		return 1;
	}


	int Valid3_1_3(uint32_t bf);
	int Valid3mm(uint32_t bf);
	void GoSubcritToMiss0(uint32_t bf, uint32_t ac);
	void GoD0F(uint32_t bf);//end all clues there

	void GoEndAll(uint32_t bf, uint32_t ac);
	void TryMiss1Subcritical();


	void GoB3Expand_1_3(uint32_t bf, uint32_t ac);
	void GoB3Expand_4_x(SP3 spe);

	void Out17(uint32_t bfb3);


};