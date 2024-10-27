#include <iostream>
#include <bitset>
#include <cmath>
#include <vector>
#include <omp.h>
#include <iomanip>
#include <random>
#include <ctime>
#include <ratio>
#include <chrono>

using namespace std;

class FP32 {
	uint32_t getsign() const noexcept {
		return data >> 31;
	}
	int32_t getexp() const noexcept {
		return (data << 1) >> 24;
	}
	uint64_t getmantissa() const noexcept {
		return data & 0x007FFFFF;
	}
	bool isfpinf() const noexcept {
		return ((getexp() == 0xFF) && getmantissa() == 0);
	}
	bool isfpnan() const noexcept {
		return ((getexp() == 0xFF) && getmantissa() != 0);
	}
	uint64_t roundDiv(uint64_t a, uint64_t qbits) const noexcept {
		uint64_t tmp = (a % (uint64_t(1) << qbits)) >= (uint64_t(1) << (qbits - 1));
		return (a >> qbits) + tmp;
	}
public:
	float example;
	uint32_t data;
	FP32() = default;
	FP32(float f) noexcept {
		data = reinterpret_cast<uint32_t*>(&f)[0];
		example = f;
	}
	FP32(const FP32& bf) noexcept {
		example = bf.example;
		data = bf.data;
	}
	FP32& operator= (const FP32& bf) noexcept {
		example = bf.example;
		data = bf.data;
		return *this;
	}
	FP32& operator= (const uint32_t& bf) noexcept {
		data = bf;
		example = *this;
		return *this;
	}
	FP32(double d) noexcept {
		example = static_cast<float>(d);
		data = reinterpret_cast<uint32_t*>(&example)[0];
	}
	operator float() const noexcept {
		uint32_t tmp = data;
		return reinterpret_cast<float*>(&tmp)[0];
	}
	FP32(uint32_t t) noexcept {
		data = t;
		example = *this;
	}
	void print() const noexcept {
		cout << bitset<32>(data) << endl;
	}

	FP32 add(const FP32 l, const FP32 r) /*noexcept*/ {
		FP32 res;
		res.example = l.example + r.example;
		uint16_t sign, normal_unit = 0;
		//NaN - checking
		int16_t el = l.getexp(), er = r.getexp();
		uint16_t ml = l.getmantissa(), mr = r.getmantissa();
		int16_t shift = el - er;
		//cout << endl << l.data << " " << r.data << " " << int(ml) << " " << int(mr) << " " << el << " " << er << endl;
		if (shift < 0) {
			res.data = er << 7;
			ml >>= -shift;
			if (7 + shift >= 0)
				normal_unit = uint16_t(1) << (7 + shift);
			//subnormal_unit = !ml;
		}
		else if (shift > 0) {
			res.data = el << 7;
			mr >>= shift;
			if (7 - shift >= 0)
				normal_unit = uint16_t(1) << (7 - shift);
			//subnormal_unit = !mr;
		}
		else {
			//overflow_unit = (er != 0);
			//res.data = (er + overflow_unit) << 7; //er << 7 if subnormals, er + 1 else
			res.data = er << 7;
			normal_unit = uint16_t(1) << 7;
			//subnormal_unit = 0;
		}
		//res.data += ((ml + mr) >> overflow_unit) + (subnormal_unit << 6);
		uint16_t mres = ml + mr + normal_unit;
		if (mres >> 7) {
			uint16_t shiftres = mres >> 7;
			res.data += (mres & 0xFF80);
			res.data += (mres & 0x007F) >> shiftres;
		}
		else {
			res.data += mres;
		}
		cout << endl << ml << " " << mr << " " << normal_unit << " " << res.data << endl;
		//res.data += ((mres & 0x7F) >> (mres >> 7)) + (mres & 0xFF80); //

		return res;
	}

	FP32 add2(const FP32 l, const FP32 r) noexcept {
		FP32 res;
		res.example = l.example + r.example;
		uint32_t sign;
		FP32 fpinf;
		FP32 fpnan;
		fpinf.data = 0x7f800000;
		fpnan.data = 0x7f800001;
		fpinf.example = INFINITY;
		fpnan.example = NAN;
		if (l.isfpnan() || r.isfpnan()) {
			return fpnan;
		}
		if (l.isfpinf() || r.isfpinf()) {
			return fpnan;
		}
		//inf - inf = nan

		int32_t el = l.getexp(), er = r.getexp();
		uint32_t ml = l.getmantissa(), mr = r.getmantissa();


		//передавать в roundDiv число - разность, а потом отбрасывать деление на числа с длиной битов меньше 1?


	}

	FP32 mul(const FP32 l, const FP32 r) noexcept {
		FP32 res;
		res.example = l.example * r.example;
		FP32 bfinf;
		FP32 bfnan;
		bfinf.data = 0x7f80;
		bfnan.data = 0x7f80 + 1;
		bfinf.example = INFINITY;
		bfnan.example = NAN;
		bool too_low = false;
		//if (l.isbfinf() || r.isbfinf()) {
		//	return bfinf;
		//}
		//if (l.isbfnan() || r.isbfnan()) {
		//	return bfnan;
		//}

		uint16_t sign = l.getsign() ^ r.getsign();
		uint16_t ml = l.getmantissa(), mr = r.getmantissa();
		int16_t el = l.getexp(), er = r.getexp();

		if ((el + er - int16_t(127)) > 0) { // >= if denormals
			res.data = (el + er - 127) << 7;
		}
		else if ((el + er - int16_t(127)) >= -7) { // matches 2^-133 // -6 -> -7
			//if ((el + er - int16_t(127)) != -7) 
			//	res.data = uint16_t(1) << (el + er - uint16_t(121));
			//else 
				res.data = 0;
			too_low = true;
		}
		else {
			res.data = 0;
			return res;
		}
		// CHECK INPUT SUBNORMALS
		uint16_t mres = ml + mr + ((ml * mr) >> 7);
		if ((mres >> 7)) {
			uint16_t shiftres = mres >> 7;
			res.data += (mres & 0xFF80);
			res.data += (mres & 0x007F) >> shiftres;
		}
		else {
			res.data += mres;
		}

		if (res.data >> 15) {
			return bfinf;
		}
		if (!res.getexp() || too_low) { //������: � ���������� ��������� ������������ �������� ��� ������� ��������� (-7) ����� ���������
			// from normal to subnormal
			//cout << endl << el + er - 127 << " " << res.data << endl;
			//if ((127 + 1 - el - er) < )
			res.data >>= (127 + 1 - el - er);
			//cout << res.data << endl;
			if (el + er - 127 + 6 >= 0) {
				res.data += uint16_t(1) << (el + er - 127 + 6);
			}
		}
		return res;
	}

	FP32 mul2(const FP32 l, const FP32 r) const noexcept { //last bit error in subnormals
		FP32 res; 
		//res.data = (l.getsign() ^ r.getsign()) << 31;
		res.data = 0;
		res.example = l.example * r.example;
		FP32 fpinf;
		FP32 fpnan;
		fpinf.data = 0x7f800000 + ((l.getsign() ^ r.getsign()) << 31);
		fpnan.data = 0x7f800001 + ((l.getsign() ^ r.getsign()) << 31);
		fpinf.example = INFINITY;
		if (l.getsign() ^ r.getsign()) fpinf.example = -INFINITY;
		fpnan.example = NAN;

		if (l.isfpinf() || r.isfpinf()) {
			return fpinf;
		}
		if (l.isfpnan() || r.isfpnan()) {
			return fpnan;
		}
		if ((l.data >> 1) == 0 || (r.data >> 1) == 0) {
			res.data = (l.getsign() ^ r.getsign()) << 31; //bad
			return data;
		}

		uint64_t mres = (l.getmantissa() + (uint64_t(l.getexp() > 0) << 23)) * (r.getmantissa() + (uint64_t(r.getexp() > 0) << 23)); //bad type (m1*m2)/2^23 = m1/2^23 * m2
		int32_t eres = l.getexp() + r.getexp() - int32_t(127) + (l.getexp() == 0 || r.getexp() == 0); // subnormal
		while ((mres < (uint64_t(1) << 46)) && (eres >= 0)) { // from subnormal to normal ZEEEEEEEE00000*****... if one arg is subnormal and m3 < 2^23
			--eres;
			mres <<= 1;
		}
		mres = roundDiv(mres, 23);
		
		
		while (mres >= (uint64_t(1) << 24)) { //instruction to count 00001***mant zeros can be used
			++eres;
			mres >>= 1;
			//mres = roundDiv(mres, 1); //works correct with simple div
		}

		if (eres > 0) {
			mres -= uint64_t(1) << 23; //as normals
			res.data += mres;
			res.data += eres << 23;
		}
		else if (eres >= -23) {
			res.data += mres/* >> -eres*/;
			res.data >>= -(eres - 1); 
			//res.data = roundDiv(res.data, -eres); //? correct? NO
		}
		else {
			res.data += (l.getsign() ^ r.getsign()) << 31; //bad
			return res;
		}
		if (res.data >> 31 || res.isfpnan()) {
			return fpinf;
		}
		res.data += (l.getsign() ^ r.getsign()) << 31;
		return res;
	}

	FP32& operator+= (const FP32& bf) noexcept {
		uint16_t sign;
		if (getsign() && bf.getsign()) sign = 1;
		else if (getsign()) {
			sign = uint16_t(1) << 15;
			*this -= bf;
			data |= sign;
			return *this;
		}
		else if (bf.getsign()) {
			*this -= bf;
			return *this;
		}
	}
	FP32& operator-= (const FP32& bf) noexcept {
		uint16_t sign;
		if (getsign() && bf.getsign()) sign = 1;
		else if (getsign()) {
			sign = uint16_t(1) << 15;
			*this -= bf;
			data |= sign;
			return *this;
		}
		else if (bf.getsign()) {
			*this -= bf;
			return *this;
		}
	}
	FP32& operator*= (const FP32& bf) noexcept {

	}
	FP32& operator/= (const FP32& bf) noexcept {

	}
	FP32 operator+ (const FP32& bf) const noexcept {
		FP32 res = *this;
		res += bf;
		return res;
	}
	FP32 operator- (const FP32& bf) const noexcept {
		FP32 res = *this;
		res -= bf;
		return res;
	}
	FP32 operator* (const FP32& bf) const noexcept {
		FP32 res = *this;
		res *= bf;
		return res;
	}
	FP32 operator/ (const FP32& bf) const noexcept {
		FP32 res = *this;
		res /= bf;
		return res;
	}
};

class Alltests {
	FP32 l;
	FP32 r;
public:
	Alltests() = default;
	bool equal_prec(uint32_t a, uint32_t b, uint32_t prec) {
		//if (a == 0 && b == 0x80000000 || a == 0x80000000 && b == 0) return true;
		for (uint32_t p = 0; p <= prec; ++p) {
			if (a - p == b || a + p == b) return true;
		}
		return false;
	}
	bool run_specific() {
		vector<uint32_t> vl = { 0xbadf8, 0x34ebd218, 0x340106e7, 0x800000, 0xbadf8 };
		vector<uint32_t> vr = { 0x40026e28, 0x801010a0, 0x0, 0x3f010130, 0x40e05790 };
		size_t from = 0;
		for (size_t i = from; i < vl.size(); ++i) {
			cout << hex << vl[i] << ", " << vr[i];
			l = vl[i];
			r = vr[i];
			//l = l.add(l, r);
			l = l.mul2(l, r);
			//l = l.add2(l, r);
			//cout << endl << float(BF16(l)) << " " << float(BF16(r)) << endl;
			//cout << l.example << " " << float(l) << endl;
			if (/*!isnan(l.example)*/ (l.example == l.example) && !equal_prec(FP32(l.example).data, l.data, 1)) {
				//cout << " ADD ERROR\n";
				cout << " MUL ERROR\n";
				cout << l.example << " expected, " << float(l) << " instead\n";
				FP32(l.example).print();
				l.print();
				return false;
			}
			cout << " PASSED\n";
		}
		return true;
	}
	bool run() {
		uint64_t lc, rc;
		for (lc = 0x00000000; lc <= 0xFFFFFFFF; lc += 7654) {
		//#pragma omp parallel for
			for (rc = 0x00000000; rc <= 0xFFFFFFFF; rc += 7654) {
				//cout << hex << lc << ", " << rc;
				//if (lc < 0x00800000 || rc < 0x00800000) continue;
				l = uint32_t(lc);
				r = uint32_t(rc);
				//l = l.add2(l, r);
				l = l.mul2(l, r);
				//cout << endl << float(BF16(l)) << " " << float(BF16(r)) << endl;
				//cout << l.example << " " << float(l) << endl;
				if (/*!isnan(l.example)*/ (l.example == l.example) && !equal_prec(FP32(l.example).data, l.data, 1)) {
					//cout << " ADD ERROR\n";
					cout << hex << endl << lc << ", " << rc;
					cout << " MUL ERROR\n";
					cout << l.example << " expected, " << float(l) << " instead\n";
					FP32(l.example).print();
					l.print();
					return false;
				}
				//cout << " PASSED\n";
			}
		}
		return true;
	}
};

template <typename T>
vector <vector <T> > mmul(const vector <vector <T> >& A, const vector <vector <T> >& B) {
	size_t m = A.size(), n = B.size(), p = B[0].size(); // A: m x n, B: n x p
	vector <vector <T> > C (m, vector<T>(p)); 

	for (size_t i = 0; i < m; ++i) {

	for (size_t j = 0; j < p; ++j) {

	for (size_t k = 0; k < n; ++k) {

	C[i][j] += A[i][k] * B[k][j];

	}
	}
	}

	return C;
}

vector <vector <double> > mmulv2(const vector <vector <double> >& A, const vector <vector <double> >& B) {
	size_t m = A.size(), n = B.size(), p = B[0].size(); // A: m x n, B: n x p
	vector <vector <double> > C(m, vector<double>(p));

	for (size_t i = 0; i < m; ++i) {
		for (size_t k = 0; k < n; ++k) {
//#pragma omp parallel for shedul static 6
			for (size_t j = 0; j < p; ++j) {


				C[i][j] += A[i][k] * B[k][j];

			}
		}
	}

	return C;
}

int main() {
	//using mat = vector <vector <double> >;
	//size_t n = 10000;
	//auto l = [](const mat& m) {
	//	for (size_t i = 0; i < m.size(); ++i) {
	//		for (size_t j = 0; j < m[i].size(); ++j) {
	//			cout << setw(10) << m[i][j];
	//		}
	//		cout << endl;
	//	}
	//};
	//auto r = [&n]() {
	//	mat m = mat(n, vector<double>(n));
	//	//srand(std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()));
	//	//srand(5);
	//	//default_random_engine generator;
	//	random_device rd;
	//	mt19937 mt(rd());
	//	uniform_real_distribution<double> dist(-100.0, 100.0);
	//	for (size_t i = 0; i < n; ++i) {
	//		for (size_t j = 0; j < n; ++j) {
	//			m[i][j] = dist(mt);
	//		}
	//	}
	//	return m;
	//};
	//
	//mat A = {
	//	{1, 5, 6},
	//	{6, 3, 8},
	//	{2, 2, 2},
	//	{1, 0, -4}
	//};
	//mat B = {
	//	{3, -8, 0, 0},
	//	{-2, -3, -3, 1},
	//	{12, 11, 6, -1}
	//};
	//
	//A = r();
	//B = r();
	//
	//cout << "ready\n";
	//auto start = chrono::high_resolution_clock::now();
	//
	//mat C = mmulv2(A, B);
	//
	//auto stop = chrono::high_resolution_clock::now();
	//auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	//cout << duration.count() << endl;
	//mat C = mmul(B, A);
	//l(C);
	bool flag = true;
	Alltests tests;
	//if (flag) flag = tests.run_specific();
	if (flag) flag = tests.run();
	cout << endl << "ENDED" << endl;
	return 0;
}