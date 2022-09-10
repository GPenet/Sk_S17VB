struct CBS {// clues band stack
	uint32_t b[3];
	uint32_t s[3];
	inline void Add(uint32_t cell) {
		b[cell / 27]++;
		s[C_stack[cell]]++;
	}
	inline int IsFilt11() {
		if (b[0] > 7 || b[1] > 6)return 1;
		if (s[0] > 7 || s[1] > 7 || s[2] > 7)return 1;
		return 0;
	}
	inline int IsFilt12() {
		if (b[0] != 6)return 1;
		if (s[0] > 6 || s[1] > 6 || s[2] > 6)return 1;
		return 0;
	}
}cbs_4, cbs_5;
struct SPB03 {// spots to first 7 clues
	BF128 v;
	uint64_t  possible_cells, all_previous_cells, active_cells;
	CBS cbs;
	uint32_t ncl;
	void Dump(uint64_t x) {
		cout << "spb03 status ncl=" << ncl << " " << x << endl;
		cout << Char54out(all_previous_cells) << " assigned" << endl;
		cout << Char54out(active_cells) << " active" << endl;
		cout << Char54out(possible_cells) << " possible" << endl;
		cout << Char64out(v.bf.u64[0]) << " 64 v" << endl;
	}
}spb_0_15[16]; 
#define UA12SIZE 3840
#define UA12BLOCS 30

struct TUASB12 {//  initial set of uas bands 1+2 
	uint64_t tua[UA12SIZE]; // 30x128
	uint32_t nua,  tdigs[UA12SIZE],ndigs[UA12SIZE];
	inline void AddInit(
		int64_t ua, uint32_t digs,  uint32_t nd) {
		tdigs[nua] = digs;
		ndigs[nua] = nd;
		tua[nua++] = ua;
	}
	void DumpInit() {
		cout << "dumpinit  tua nua=" << nua << endl;
		for (uint32_t iua = 0; iua < nua; iua++) {
			cout << Char2Xout(tua[iua]) << " i="
				<< iua << " " << (tua[iua] >> 59)
				<< " " <<_popcnt64(tua[iua] & BIT_SET_2X);
			cout << " digs " << Char9out(tdigs[iua]) << " " << ndigs[iua] << endl;;
		}
	}
	void DumpShort(const char * lib) {
		cout << lib  << "tuasb12 tua nua=" << nua << endl;
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
		void Dumpv(int lim = 32) {
			cout << Char32out(v0.bf.u32[0]) << " v0" << endl;
			for (int i = 0; i < 54; i++)
				cout << Char32out(vc[i].bf.u32[0]) << " cell=" << i << endl;
		}
		void Dump(int lim = 128) {
			for (int i = 0; i < lim; i++)
				if (v0.On(i))
					cout << Char54out(t[i]) << " " << i <<" "<<_popcnt64(t[i]) << endl;
				else return;
		}
	};

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
	int Build_tb128();

#define UACNBLOCS 15
	// status after 3 clues (C)
	TUVECT tc128[UACNBLOCS];// max start 10*128=1280 uas 
	uint32_t nc128, ncblocs, ntc128[UACNBLOCS];
	void InitC() {
		memset(ntc128, 0, sizeof ntc128);
		nc128 = ncblocs = 0;
		for (int i = 0; i < UACNBLOCS; i++)tc128[i].Init();
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
		register uint32_t cell;
		while (bitscanforward64(cell, R)) {
			R ^= (uint64_t)1 << cell; //clear bit
			myvc[cell].clearBit(ir);
		}
	}
	int Build_tc128();
	inline int IsNotRedundant(uint64_t u) {
		register uint64_t nu = ~u;
		for (uint32_t i = 0; i < ntc128[0]; i++)
			if (!(tc128[0].t[i] & nu)) return 0;
		return 1;
	}

}t54b12;

struct T54G2 {//
	struct G2VECT {//  128 guas2 and vectors  
		BF128 v0, vc[54], v81[81];
		//vectors 128 {base, cells, i81s}
		uint64_t t[128];// ua54 + i81
		void Init() {
			v0.SetAll_0();
			memset(vc, 255, sizeof vc);
			memset(v81, 255, sizeof v81);
		}
	};
#define NG2BLOCS6 30
	G2VECT t128[NG2BLOCS6];// max start 20*128=2560 uas 
	uint32_t n128, nblocs, nt128[NG2BLOCS6];
	void Init() {
		memset(nt128, 0, sizeof nt128);
		n128 = nblocs = 0;
		for (int i = 0; i < NG2BLOCS6; i++)t128[i].Init();
	}
	inline void Add(uint64_t u) {
		if (n128 >= NG2BLOCS6 * 128) return;
		register uint32_t bloc = n128 >> 7, ir = n128 - (bloc << 7);
		n128++; nblocs = bloc; nt128[bloc]++;
		register G2VECT& myb = t128[bloc];
		myb.v0.setBit(ir);
		BF128* myvc = myb.vc,*myvi=myb.v81;
		register uint64_t R = u;
		myb.t[ir] = R;
		myvi[R>>56].clearBit(ir);
		R &= BIT_SET_54;// clear i81
		uint32_t cell;
		while (bitscanforward64(cell, R)) {
			R ^= (uint64_t)1 << cell; //clear bit
			myvc[cell].clearBit(ir);
		}
	}
	void BuildG2();

}t54g2;
struct T54G3 {//
	struct G3VECT {//  128 guas2 and vectors 
		BF128 v0, vc[54], v81[81];
		uint64_t t[128];// ua54 + i81
		void Init() {
			v0.SetAll_0();
			memset(vc, 255, sizeof vc);
			memset(v81, 255, sizeof v81);
		}
	};

#define NG3BLOCS 15
	// initial status after harvest plus fresh guas
	G3VECT t128[NG3BLOCS];// max start 15*128=1920 uas 
	uint32_t n128, nblocs, nt128[NG3BLOCS];
	void Init() {	
		memset(nt128, 0, sizeof nt128);
		n128 = nblocs = 0;
		for (int i = 0; i < NG3BLOCS; i++)t128[i].Init();
	}
	inline void Add(uint64_t u) {
		if (n128 >= NG3BLOCS * 128) return;
		register uint32_t bloc = n128 >> 7, ir = n128 - (bloc << 7);
		n128++; nblocs = bloc; nt128[bloc]++;
		register G3VECT& myb = t128[bloc];
		myb.v0.setBit(ir);
		BF128* myvc = myb.vc, * myvi = myb.v81;
		register uint64_t R = u;
		myb.t[ir] = R;
		myvi[R >> 56].clearBit(ir);
		R &= BIT_SET_54;// clear i81
		uint32_t cell;
		while (bitscanforward64(cell, R)) {
			R ^= (uint64_t)1 << cell; //clear bit
			myvc[cell].clearBit(ir);
		}
	}
	void BuildG3();

}t54g3;




#define NG2_128 6
#define NG23_128 10
struct G_128 {//vector 256 chunk level gua2 gua3
	BF128 v, vc[54], vp[27];
	uint32_t bits[128];//bit  i27 or i9
	uint32_t pat27[128];// i27 or i9

	uint32_t nv;
	void Init() {
		nv = 0;
		memset(&v, 0, sizeof v);
		memset(vc, 255, sizeof(vc));
		memset(vp, 255, sizeof(vp));
	}

	void Add(BF128 gx) { // gx is a still valid gua2 gua3 
		if (nv >= 128)return; // should never be
		int ix = 0, nx = nv;
		v.Set(nx);
		bits[nv] = gx.bf.u32[2];
		register uint64_t V = gx.bf.u64[0];
		register uint32_t x;
		while (bitscanforward64(x, V)) {
			V ^= (uint64_t)1 << x;// clear bit
			vc[x].clearBit(nx);
		}
		register uint32_t i27 = gx.bf.u32[3];
		pat27[nv] = i27;
		if (i27 < 27)vp[i27].clearBit(nx);
		nv++;
	}
	void GetActive0(BF128 vw, uint32_t& gbf,
		uint32_t* tg, uint32_t& ntg) {
		for (register uint32_t i = 0; i < ntg; i++)
			vw &= vp[tg[i]];//apply olds
		register uint64_t V = vw.bf.u64[0];
		register uint32_t  x;
		while (bitscanforward64(x, V)) {
			V ^= (uint64_t)1 << x;
			gbf |= bits[x];
			register uint32_t i27 = pat27[x];
			tg[ntg++] = i27;
			vw &= vp[i27];
			V &= vp[i27].bf.u64[0];
		}
		V = vw.bf.u64[1];
		while (bitscanforward64(x, V)) {
			V ^= (uint64_t)1 << x;
			gbf |= bits[x + 64];
			register uint32_t i27 = pat27[x + 64];
			tg[ntg++] = i27;
			V &= vp[i27].bf.u64[1];
		}
	}

	uint32_t Apply(uint32_t* tc, uint32_t ntc) {
		if (!nv) return 0;
		BF128 vw=v;
		uint32_t g2 = 0;
		for (uint32_t i = 0; i < ntc; i++) {// apply
			uint32_t cell = tc[i];
			vw &= vc[cell];
		}
		for (uint32_t i = 0; i < 2; i++) {// extract
			register uint64_t V = vw.bf.u64[i];
			register uint32_t* b = &bits[64 * i], x;
			while (bitscanforward64(x, V)) {
				V ^= (uint64_t)1 << x;
				g2 |= b[x];
			}
		}
		return g2;
	}
	void ApplyMore(uint32_t* tc, uint32_t ntc,
		uint32_t* tmore, uint32_t& ntmore, uint32_t asb3) {
		if (!nv) return;
		register uint32_t Rasb3 = asb3;
		BF128 vw=v;
		for (uint32_t i = 0; i < ntc; i++) {// apply
			uint32_t cell = tc[i];
			vw &= vc[cell];
		}
		register uint64_t* pV = vw.bf.u64;
		for (uint32_t i = 0; i < 4; i++) {// extract
			register uint64_t V = pV[i];
			register uint32_t* b = &bits[64 * i], x;
			while (bitscanforward64(x, V)) {
				V ^= (uint64_t)1 << x;
				register uint32_t U = b[x];
				if (!(Rasb3 & U))	tmore[ntmore++] = U;
			}
		}
	}
	void Dump() {
		for (uint32_t i = 0; i < nv; i++) {
			cout << Char27out(bits[i]) << " ";
			for (int j = 0; j < 54; j++)
				if (vc[j].On(i)) cout << ".";
				else cout << "1";
		}
	}

}gx_128x[NG23_128]; //g2_128[4], g3_128[2], gm_256[2];
struct G2_64Add {
	uint64_t u[64];//ua
	uint32_t p[64];// i27 or i9 or pat
	uint32_t nv;
	void Add(BF128 w) {
		if (nv >= 64)return;
		u[nv] = w.bf.u64[0];
		p[nv++] = w.bf.u32[2];
	}
	void Getactive(uint64_t bf, uint32_t& gbf,
		uint32_t* tg, uint32_t& ntg) {
		for (uint32_t i = 0; i < nv; i++) {
			if (u[i] & bf) continue;
			register uint32_t P = p[i], bit = 1 << P;
			if (gbf & bit) continue;
			gbf |= bit;
			tg[ntg++] = P;
		}
	}
	void Dump() {
		cout << "dump G2_64Add" << endl;
		for (uint32_t i = 0; i < nv; i++) {
			cout << Char27out(p[i]);
			cout << " " << Char54out(u[i]) << endl;

		}
	}
}g2_64, g3_64, gm_64;
struct G_128Handler {
	G_128* g2_128, * g3_128;// , * gm_128;
	uint32_t n2,n3, nm;
	BF128 vcl[7][NG23_128];//6 clues  
	uint64_t bf12;
	uint32_t g2, g3, tg2[27], tg3[9], ntg2, ntg3;
	G_128Handler() {
		g2_128 = gx_128x; 
		g3_128 = &gx_128x[NG2_128];
		//gm_128 = &g3_128[2];
	}
	inline void GetStart(uint32_t n2e, uint32_t n3e, uint32_t nme) {
		n2 = n2e;  n3 = n3e; nm = nme;
		BF128* V = vcl[0];
		for (int i = 0; i < NG23_128; i++) V[i] = gx_128x[i].v;
		g2_64.nv = 0; g3_64.nv = 0; gm_64.nv = 0;
	}

	inline void Add128x(BF128 w, uint32_t cc, int add = 0) {
		if (add) {
			if (cc == 2) g2_64.Add(w);
			else if (cc == 3) g3_64.Add(w);
			else gm_64.Add(w);
		}
		if (cc < 4) {
			w.bf.u32[2] = 1 << w.bf.u32[2];// switch to bit field
			if (cc == 2) {
				//cout << "add 2 in g2_256 n2=" << n2 << endl;
				//if (n2++ < 256)g2_256[0].Add(w);
				//else g2_256[1].Add(w);
			}
			//else { g3_256.Add(w); }
			return;
		}
		//int i = nm++ >> 7;
		//if (nm++ < 256)gm_128[i].Add(w);
	}
/*
	void ApplyMore(uint32_t* tc, uint32_t ntc,
		uint32_t* tmore, uint32_t& ntmore, uint32_t asb3) {
		int nb = nm++ >> 7;
		for(int i=0;i<nb;i++)
			gm_128[i].ApplyMore(tc, ntc, tmore, ntmore, asb3);
		if (gm_64.nv) {
			for (uint32_t i = 0; i < gm_64.nv; i++) {
				if (gm_64.u[i] & bf12) continue;
				register uint32_t P = gm_64.p[i];
				if (asb3 & P) continue;
				tmore[ntmore++] = P;
			}
		}
	}*/


	inline void NewVcl(int ncl, int clue) {
		register int i1 = ncl - 7;
		register BF128* V1 = vcl[i1], * V2 = vcl[i1 + 1];
		for(int i=0;i< NG23_128;i++)
			V2[i] = V1[i] & gx_128x[i].vc[clue];
	}
	void GetActive(int nclues, uint64_t bf) {
		BF128* V2 = vcl[nclues - 6];
		bf12 = bf;
		g2 = g3 = ntg2 = ntg3 = 0;
		for (int i = 0; i < NG2_128; i++)		if (V2[i].isNotEmpty())
			gx_128x[i].GetActive0(V2[i], g2, tg2, ntg2);
		for (int i = NG2_128; i < NG23_128; i++)		if (V2[i].isNotEmpty())
			gx_128x[i].GetActive0(V2[i], g3, tg3, ntg3);
		if (g2_64.nv)g2_64.Getactive(bf, g2, tg2, ntg2);
		if (g3_64.nv)g3_64.Getactive(bf, g3, tg3, ntg3);
	}
	void Dumpvcl(int nclues) {
		cout << "Dump vcl nclues=" << nclues << endl;
		BF128* V2 = vcl[nclues - 6];
		uint64_t* v64 = V2[0].bf.u64;
		for (int i = 0; i < NG23_128; i++) {
			register uint64_t v = v64[i];
			if (v)	cout << Char64out(v) << " i=" << i
				<< " cpt=" << _popcnt64(v) << endl;

		}
	}

}g_128h;

struct CHUNK3B {// storing 64 uas and vectors 3 bands
	uint64_t tu54[64], v0,
		vc[54], //cell for tu54
		vclean, // v0 after clean group
		*vbx,//valid for b3 x
		nt;
	uint32_t tu27[64];// index 0-80 (g2 g3) or  pattern
	inline void Init() {
		nt = 0, v0 = vclean = 0;
		memset(vc, 255, sizeof vc);
	}
	inline void Add(uint64_t u54, uint32_t u27) {//add a new ua
		if (nt >= 64) return;// safety should never be
		//if (1 && u27 == 1) {
		//	cout <<Char54out(u54)<< " add nt=" << nt << endl;
		//}
		uint64_t bit = (uint64_t)1 << nt;
		tu27[nt] = u27;
		tu54[nt++] = u54;
		v0 |= bit;
		uint32_t cc54;// build cells vectors
		register  uint64_t Rw = u54;
		while (bitscanforward64(cc54, Rw)) {
			Rw ^= (uint64_t)1 << cc54;// clear bit
			vc[cc54] ^= bit;
		}
	}
	inline void ApplyClean(uint32_t* tc, uint32_t ntc) {
		register uint64_t V = v0;
		for (uint32_t i = 0; i < ntc; i++)
			V &= vc[tc[i]];
		vclean = V;
	}

	void Get2x(BF128* td, uint32_t& ntd,
		uint32_t* tc, uint32_t ntc) {
		register uint64_t V = v0;
		uint32_t ir;
		for (uint32_t i = 0; i < ntc; i++)
			V &= vc[tc[i]];
		vclean = V;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			BF128& w = td[ntd++];
			w.bf.u64[0] = tu54[ir];
			w.bf.u32[2] = tu27[ir];
		}
	}

	void GetAddClean(uint64_t* td54, uint32_t* td27, uint32_t& ntd) {
		ntd = 0;
		register uint64_t V = vclean;
		register uint32_t ir;
		while (bitscanforward64(ir, V)) {
			V ^= (uint64_t)1 << ir;
			td27[ntd] = tu27[ir];
			td54[ntd++] = tu54[ir];
		}
	}
	void DebugC2(int all = 1) {
		for (uint32_t i = 0; i < nt; i++) {
			cout << Char54out(tu54[i]) << " " << tu27[i] << "\ti=" << i << endl;
		}
		if (!all) return;
		cout << Char64out(v0) << " v0" << endl;
		for (int i = 0; i < 54; i++)
			if ((vc[i] & v0) != v0)
				cout << Char64out(vc[i] & v0) << " cell=" << i << endl;
	}
	void DebugC2Clean() {
		for (uint32_t i = 0; i < nt; i++) {
			uint64_t bit = (uint64_t)1 << i;
			if (vclean & bit)
				cout << Char54out(tu54[i]) << " " << tu27[i] << " i=" << i << endl;
		}
	}
	void DebugMore(int all = 1) {
		for (uint32_t i = 0; i < nt; i++) {
			cout << Char54out(tu54[i]) << " ";
			cout << Char27out(tu27[i]) << " i=" << i << endl;
		}
		if (!all) return;
		cout << Char64out(v0) << " v0" << endl;
		for (int i = 0; i < 54; i++)
			if ((vc[i] & v0) != v0)
				cout << Char64out(vc[i] & v0) << " cell=" << i << endl;
	}
};
struct CHUNK1B {//storing 64 uas and vectors band 3s
	uint64_t v0, vc[27], nt;
	uint32_t tua[64];
	inline void Init() {
		nt = 0, v0 = 0;
		memset(vc, 255, sizeof vc);
	}
	inline void Add(uint32_t v) {//add a new ua
		if (nt >= 64) return;// safety should never be

		{//__ check redundancy
			register uint32_t vn = ~v;
			for (uint64_t i = 0; i < nt; i++)
				if (!(tua[i] & vn))return; // == or subset
		}

		uint64_t bit = (uint64_t)1 << nt;
		v &= BIT_SET_27;//no extra bit
		tua[nt++] = v;
		v0 |= bit;
		uint32_t cc27;// build cells vectors
		register  uint32_t Rw = v;
		while (bitscanforward(cc27, Rw)) {
			Rw ^= 1 << cc27;// clear bit
			vc[cc27] ^= bit;
		}
	}


	inline uint64_t ApplyXY(uint32_t* tcells, uint32_t ntcells) {
		if (!nt) return 0;
		uint64_t w = v0;
		for (uint32_t i = 0; i < ntcells; i++)
			w &= vc[tcells[i]];
		return w;
	}
	void Debug(int all = 1) {
		for (uint32_t i = 0; i < nt; i++) {
			cout << Char27out(tua[i]) << " i=" << i << endl;
		}
		if (!all) return;
		cout << Char64out(v0) << " v0" << endl;
		for (int i = 0; i < 27; i++)
			if ((vc[i] & v0) != v0)
				cout << Char64out(vc[i] & v0) << " cell=" << i << endl;
	}
}b3direct;
struct CHUNKS_HANDLER {
	BF128 	t2digs[200]; //  has band 3 after 2 digits

#define CSIZE 99
	CHUNK3B c2[CSIZE + 1], c3[CSIZE + 1],
		c4[CSIZE + 1], c5[CSIZE + 1], cmore[CSIZE + 1];
	CHUNK1B band3[2];// up to 81 band3 up to 128 others
	uint32_t ic2, ic3, ic4, ic5, icmore, iband3;
	//___________ valid to process t4 5 6 more
	uint32_t t12[1100], nt12,nt2digs;
	void Init() {
		ic2 = ic3 = ic4 = ic5 = icmore = iband3 = 0;
		c2[0].Init();  c3[0].Init();
		c4[0].Init();  c5[0].Init(); cmore[0].Init();
		band3[0].Init(); band3[1].Init();
	}
	inline int GetC2Count() { return 64 * ic2 + (int)c2[ic2].nt; }
	inline void Add128(BF128 w54, uint32_t cc) {
		switch (cc) {
		case 2:
			//if(Check2(w54 ))return;
			if (ic2 == CSIZE && c2[CSIZE].nt > 63) return;
			if (c2[ic2].nt > 63)c2[++ic2].Init();
			c2[ic2].Add(w54.bf.u64[0], w54.bf.u32[2]);
			return;
		case 3:
			//if (Check3(w54))return;
			if (ic3 == CSIZE && c3[CSIZE].nt > 63) return;
			if (c3[ic3].nt > 63)c3[++ic3].Init();
			c3[ic3].Add(w54.bf.u64[0], w54.bf.u32[2]);
			return;
		case 4:
			//if (Check4(w54))return;
			if (ic4 == CSIZE && c4[CSIZE].nt > 63) return;
			if (c4[ic4].nt > 63)c4[++ic4].Init();
			c4[ic4].Add(w54.bf.u64[0], w54.bf.u32[2]);
			return;
		case 5:
			if (ic5 == CSIZE && c5[CSIZE].nt > 63) return;
			if (c5[ic5].nt > 63)c5[++ic5].Init();
			c5[ic5].Add(w54.bf.u64[0], w54.bf.u32[2]);
			return;
		default:
			if (icmore == CSIZE && cmore[CSIZE].nt > 63) return;
			if (cmore[icmore].nt > 63)cmore[++icmore].Init();
			cmore[icmore].Add(w54.bf.u64[0], w54.bf.u32[2]);
			return;
		}

	}

	void Addc2(uint64_t ua12, uint32_t i81) {
		uint64_t ua54 = (ua12 & BIT_SET_27) |
			((ua12 & BIT_SET_B2) >> 5);
		if (ic2 == CSIZE && c2[CSIZE].nt > 63) return;
		if (c2[ic2].nt > 63)c2[++ic2].Init();
		c2[ic2].Add(ua54, i81);

	}
	void Addband3(uint32_t w) {
		if (iband3 == 1 && band3[1].nt > 63) return;
		if (band3[iband3].nt > 63)band3[++iband3].Init();
		band3[iband3].Add(w);
	}

	inline void ApplyClean(uint32_t* tc, uint32_t ntc) {
		for (uint32_t i = 0; i <= ic2; i++)
			c2[i].ApplyClean(tc, ntc);
		for (uint32_t i = 0; i <= ic3; i++)
			c3[i].ApplyClean(tc, ntc);
		for (uint32_t i = 0; i <= ic4; i++)
			c4[i].ApplyClean(tc, ntc);
		for (uint32_t i = 0; i <= ic5; i++)
			c5[i].ApplyClean(tc, ntc);
		for (uint32_t i = 0; i <= icmore; i++)
			cmore[i].ApplyClean(tc, ntc);
	}

	void DebugAll(int full = 0) {
		cout << "chunkh debug all"
			<< "\nic2/nc2 " << ic2 << " " << c2[ic2].nt
			<< "\nic3/nc3 " << ic3 << " " << c3[ic3].nt
			<< "\nic4/nc4 " << ic4 << " " << c4[ic4].nt
			<< "\nic5/nc5 " << ic5 << " " << c5[ic5].nt
			<< "\nicm/ncm " << icmore << " " << cmore[icmore].nt << endl;
		if (!full) return;
		if (full == 2) {
			cout << "band 3" << endl;
			for (uint32_t i = 0; i <= iband3; i++)
				band3[i].Debug(0);
			cout << "guas 2 cells" << endl;
			for (uint32_t i = 0; i <= ic2; i++)
				c2[i].DebugC2(0);
		}
		cout << "guas 3 cells" << endl;
		for (uint32_t i = 0; i <= ic3; i++)
			c3[i].DebugC2(0);
		cout << "guas more cells" << endl;
		for (uint32_t i = 0; i <= ic4; i++)
			c4[i].DebugMore(0);
		for (uint32_t i = 0; i <= ic5; i++)
			c5[i].DebugMore(0);
		for (uint32_t i = 0; i <= icmore; i++)
			cmore[i].DebugMore(0);

	}

	void C2Status() {
		for (uint32_t i = 0; i <= ic2; i++)
			c2[i].DebugC2(0);
	}
	void C4Status() {
		for (uint32_t i = 0; i <= ic4; i++)
			c4[i].DebugMore(0);
	}

	void C3StatusClean() {
		cout << "status for c3 clean" << endl;
		for (uint32_t i = 0; i <= ic3; i++)
			c3[i].DebugC2Clean();
	}
	void Debug27() {
		cout << "chunk extract status n=" << nt12 << endl;
		for (uint32_t i = 0; i < nt12; i++)
			cout << Char27out(t12[i]) << endl;
	}
	void Stats() {
		cout << " chunks status ";
		cout << " ic2 " << ic2 << ";" << c2[ic2].nt;
		cout << " ic3 " << ic3 << ";" << c3[ic3].nt;
		cout << " ic4 " << ic4 << ";" << c4[ic4].nt;
		cout << " ic5 " << ic5 << ";" << c5[ic5].nt;
		cout << " icmore " << icmore << ";" << cmore[icmore].nt;
		cout << endl;
	}
}chunkh;

struct MINCOUNT {
	uint32_t mini_bf1, mini_bf2, mini_bf3, mini_triplet,
		critbf, pairs27;
	uint32_t all_used_minis, mincount, minplus;
	inline void SetMincount() {// after direct setting minis
		all_used_minis = mini_bf1 | mini_triplet;
		mini_triplet &= ~mini_bf1;// count only triplets with no pair
		mincount = _popcnt32(all_used_minis) + _popcnt32(mini_bf3);
		mini_bf1 &= ~mini_bf2;// now pure one pair
		mini_bf2 &= ~mini_bf3;// now pure 2 pairs

		// set up pair + triplet bitfield
		if (mini_triplet) {// must add triplet minirow
			for (int i = 0, bit = 1, field = 7; i < 9; i++, bit <<= 1, field <<= 3)
				if (mini_triplet&bit)
					critbf |= field;
		}
		minplus = mincount;
	}
	GINT64_t  Count_per_stack() {
		GINT64 cc; cc.u64 = 0;
		for (int i = 0, st = 0111; i < 3; i++, st <<= 1) {
			cc.u16[i] = _popcnt32(all_used_minis&st) +
				_popcnt32(mini_bf3&st);
		}
		return cc;
	}
	void Status(const char * lib) {
		cout << lib << "critical Status mincount =" << mincount << " minplus=" << minplus << endl;
		cout << Char27out(critbf) << " critical bf" << endl;
		cout << Char27out(pairs27) << " pairs 27" << endl;
		cout << Char9out(mini_bf1) << "     minis bf1" << endl;
		cout << Char9out(mini_bf2) << "     minis bf2" << endl;
		cout << Char9out(mini_bf3) << "     minis bf3" << endl;
		cout << Char9out(mini_triplet) << " mini triplets" << endl << endl;

	}
};


// standard first band (or base any band band)
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
	void MorphUas()	;
	void InitC10(int i);// known mode


	void PrintStatus();
};
struct STD_B3 :STD_B416 {// data specific to bands 3
	BF128 tua128[1000], tua128_b2[1000];
	struct G {
		BF128 gsocket2, gsocket3;// active i81 mode 81 bits
		int pat2[81], pat3[81]; // storing ua bitfields
		int ua2_imini[81], ua3_imini[81],	ua2pair27[81];
	}g;


	uint32_t ntua128,ntua128_b2;
	int minirows_bf[9];
	int triplet_perms[9][2][3];
	uint32_t i_27_to_81[27], i_9_to_81[9]; //band i81 for guas guas3
	//_______________ handling mincount
	MINCOUNT smin,sminr;
	GINT64  stack_count;

	//_______________________

	void InitBand3(int i16, char * ze, BANDMINLEX::PERM & p);

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
	void Pat_to_digits(int bf,int* tcells, int* tdigits, int& nt) {
		register int cell;
		nt = 0;
		while (bitscanforward(cell, bf)) {
			bf ^= 1 << cell;
			tcells[nt] = cell;
			tdigits[nt++] = band0[cell];
		}
	}
	int Is_Pat_For_Me(int* tcells, int* tdigits, int nt) {
		for (int i = 0; i < nt; i++) {
			if (band0[tcells[nt]] != tdigits[nt]) return 0;
		}
		return 1;
	}
	/*
	inline void Insert2(uint32_t i81) {
		//tua2_3[ntua2_3++] = guas.ua_pair[i81];
		register uint32_t	bit = guas.ua2bit[i81];
		smin.mini_bf3 |= smin.mini_bf2&bit;
		smin.mini_bf2 |= smin.mini_bf1&bit;
		smin.mini_bf1 |= bit;
		smin.critbf |= guas.ua_pair[i81];
		smin.pairs27 |= guas.ua2pair27[i81];
	}
	inline void Insert3(uint32_t i81) {
		//tua2_3[ntua2_3++] = guas.ua_triplet[i81];
		smin.mini_triplet |= guas.ua3bit[i81];
	}	*/

	void PrintB3Status();
};

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
	//uint32_t nua;// nua_start, nua_end;

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
	//uint32_t nua;// nua_start, nua_end, nua;
}tsgua3[81];

struct GEN_BANDES_12 {// encapsulating global data 
	STD_B3 bands3[512];
	int modeb12, go_back, diagmore, diagbug, ip20,
		it16, it16_2, imin16_1, imin16_2, imin16_3;
	int i1t16, i2t16, i3t16, maxnb3; // index 416 ordered in increasing size of valid clues 6
	char zsol[82], rband2[28];
	int grid0[81], tc[6], ntc;
	int gcheck[82], ib2check, ib3check;
	int skip, last;// restart point; last entry in the batch
	uint64_t   nb12;
	BANDMINLEX::PERM t_auto_b1[108], // maxi is 107excluding start position
		t_auto_b1b2[108], t_auto_b2b1[108],
		pband2, pband3, pcheck2, pcheck3;
	int n_auto_b1, n_auto_b1b2, n_auto_b2b1;
	int cold[9], coldf[9], rowd[6], boxd[6], rowdb3[3], boxdb3[3]; //free digits 
	//_________________ gangster 
	int gangcols[9];// 9 bits fields band3 for each colum (see valid band2)
	//int gang[9][3]; // gangcols expanded (buildgang ) 3 digits
	int gangb12[9]; // digit bf bands 12 per column
	//int   *gang27; // redefines gang[9][3] as 27 integer
	//int   gang_digits_cols[9][3];// active cols for a given digit
	//____________structs hosting the 81 GUA entries


	int GET_I81_G2(int digs, uint32_t pat) {
		register  uint32_t A = pat,
			B = (A | (A >> 9) | (A >> 18)) & 0777; // all columns
		uint32_t dstack;
		bitscanforward(dstack, B);
		dstack = (dstack / 3) * 27;// now stack0 0, 27 57 chunks of 27 i81
		for (uint32_t i = dstack; i < dstack + 27; i++)
			if (digs == tsgua2[i].digs) return i;
		return 0; // would be a bug
	}
	int GET_I81_G2_4(int digs, uint32_t pat) {
		register  uint32_t A = pat,
			B = (A | (A >> 9) | (A >> 18)) & 0777; // all columns
		// keep the stack with 2 columns
		for (int i = 0, mask = 7; i < 3; i++, mask <<= 3) {
			if (_popcnt32(B&mask) == 2) {
				B = mask;
				break;
			}
		}
		uint32_t dstack; // now same as pat 2 cells
		bitscanforward(dstack, B);
		dstack = (dstack / 3) * 27;// now stack0 0, 27 57 chunks of 27 i81
		for (uint32_t i = dstack; i < dstack + 27; i++)
			if (digs == tsgua2[i].digs) return i;
		return 0; // would be a bug
	}

	int GET_I81_G3(int digs, uint32_t pat) {
		register  uint32_t A = pat,
			B = (A | (A >> 9) | (A >> 18)) & 0777; // all columns
		uint32_t dstack;
		bitscanforward(dstack, B);
		dstack = (dstack / 3) * 27;// now stack0 0, 27 57 chunks of 27 i81
		for (uint32_t i = dstack; i < dstack + 27; i++)
			if (digs == tsgua3[i].digs) return i;
		return 0; // would be a bug
	}

	// __________________________  primary UAs tables and creation of such tables
	uint64_t  *ptua2;// pointer to current table cycle search 2/3
	uint32_t   nua2; // nua2 for cycle search 2/3  
	//================== bands 3 and gangster band 3 analysis
	int nband3;
	int   tcolok[2], ncolok;

	int ngua6_7, c1, c2, band, floors, digp, i81;
	uint64_t wua0, ua;// partial gua ua to check
	uint64_t tuacheck[100], tua_for_check[500];
	uint32_t uadigs_for_check[500], nua_for_check, nua_check;
	//================ A creating a catalogue for the 17 search 
	//sorted increasing number of valid bands 6 clues

	GEN_BANDES_12() {
		//gang27 = gang[0];
		//InitialSockets2Setup();
		//InitialSockets3Setup();
	}
	//void InitialSockets2Setup();// batch level
	//void InitialSockets3Setup();// batch level
	//================================= functions
	void GetStartB2(int i); // one of the 20 starts 
	void Start(int mode = 0);
	void NewBand1(int iw);
	int Band2Check();
	int Band3Check();
	void Find_band2B();
	int ValidBand2();
	void ValidInitGang();
	void Find_band3B(int m10 = 1);
	void Find_band3B_pass1(int m10=1);
	//============= loops control for UAs 5;6;7 digits collection (collect more=
	int iband, ibox, iminirow, ibox2, iminirow2, pat1, pat2, ncells;
	int tcells[6], tcols[6];
	int bcols[2][9], mycols[9], myfloors;
	uint64_t mybf;
};

struct G17B {// hosting the search in 6 6 5 mode combining bands solutions
	G17B();// initial tasks all commands

	BF128 p17diag;// known 17 pattern for tests
	int mincluesb3,
		b3lim,		 aigstop, aigstopxy,
		//iretb1,doloopnotok,
		npuz, a_17_found_here;
	int  debug17, debug17_check, diag, diagbug, debugb3,
		is_test_on,ng2,ng3;
	//uint32_t	iband1,iband2, step1count;
	// band in use
	int grid0[81];

	//____gangsters, brute force,sockets setup
	//_________________ gangster 
	//int gangcols[9];// 9 bits fields band3 for each colum (see valid band2)
	int gang[9][3]; // gangcols expanded (buildgang ) 3 digits
	//int gangb12[9]; // digit bf bands 12 per column
	int* gang27; // redefines gang[9][3] as 27 integer
	int   gang_digits_cols[9][3];// active cols for a given digit

	//______sockets common to  all bands 3  
	BF128 gsock2, gsock3;

	//====== data for band expansion
	uint32_t nexp, bnua;
	uint64_t btua[300],start_active, b1cpt[8], b2cpt[8],b1cptdiag;
	uint32_t nmybi2t, nmyokt,nbi2_1,nbi2_2,
		nvb1,nvb1steps[10],
		nzs1_6,nzs2_6, nzs1_5, nzs2_5;
	//======= status after step 2 in band 2 then 2 uas in band 1
	uint64_t tusb1[2000], tusb2_12[1000], tusb2[1000], temptyb1[500];
	uint32_t ntusb1, ntusb2, ntusb2_12, nemptyb1;
	uint64_t  temptyb2[500];
	uint32_t  ntua_128, ntua_256, nemptyb2;
	uint64_t fb1,acb1, fb2, fb12, acb2a,  acb2, acb12;



	//============================ b12 no more uas to test
	uint64_t wb12bf, wb12active,myua;
	GINT64  stack_count_step, stack_count, stack_countf;
	uint32_t tclues[40],tb3[256],*tcluesxy; 
	int nclues_step, nclues,ntb3;
	BF128 bands_active_pairs, bands_active_triplets,
		valid_vect;
	//============  go band3
	STD_B3* myband3;
	uint32_t fstk, andmiss1, noutmiss1, wactive0;
	uint32_t free1, free2, free3;
	uint32_t tua128_b3[1000], ntua128_b3;
	BF128 wg46;
	uint32_t cur_ib;
	uint32_t tcluesb12[20], ncluesb3x;
	uint32_t   nmiss;
	uint32_t uasb3_1[2000], uasb3_2[2000], uas_in[2000],
		nuasb3_1, nuasb3_2, nuas_in, b3_andout;
	MINCOUNT smin;

	
	inline void AddB3_2_to_B3_1() {
		memcpy(&uasb3_1[nuasb3_1], uasb3_2, nuasb3_2 * 4);
		nuasb3_1 += nuasb3_2;
	}
	
	
	
	
	//=====================process for a new band 2 / set of bands 3
	void Start();// standard entry
	void StartPrint();// standard entry
	void StartKnown();// entry for a known 17
	void StartInit();// initial task gangster set up 
	void StartInitDebug();// initial task gangster set up 
	void UaCollector();
	inline void Adduab12(uint32_t digs, uint32_t nd);
	void FirstUasCollect();
	void SecondUasCollect();
	void UasCollect4box();
	void UasCollect6_7();

	void StartAfterUasHarvest();
	void SAH_BuildGuas2();
	void SAH_BuildGuas3();

	inline int BuildGua(BF128& w);
	inline void BuildGua(BF128& w, int cc);
	void Guas2Collect();
	void Guas2CollectG2();
	void Guas2CollectG3();
	void Guas2CollectG3_4d();
	void Guas2CollectG3_5d();
	void Expand_03();

	int SetupExpand_46();
	void Expand_46();

	int SetupExpand_7p();
	void Init7p_guas();

	int nt4ok, okcheck;// for known
	// bands 1+2 valid epansion


};