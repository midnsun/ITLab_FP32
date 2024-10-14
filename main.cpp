#include <iostream>
#include <bitset>
#include <cmath>

using namespace std;

class BF16 {
	uint16_t getsign() const noexcept {
		return data >> 15;
	}
	int16_t getexp() const noexcept {
		return (data << 1) >> 8;
	}
	uint16_t getmantissa() const noexcept {
		return data & 0x7F;
	}
	bool isbfinf() const noexcept {
		return ((getexp() == 0xFF) && getmantissa() == 0);
	}
	bool isbfnan() const noexcept {
		return ((getexp() == 0xFF) && getmantissa() != 0);
	}
public:
	float example;
	uint16_t data;
	BF16() = default;
	BF16(float f) noexcept {
		data = reinterpret_cast<uint16_t*>(&f)[1];
		example = f;
	}
	BF16(const BF16& bf) noexcept {
		example = bf.example;
		data = bf.data;
	}
	BF16& operator= (const BF16& bf) noexcept {
		example = bf.example;
		data = bf.data;
		return *this;
	}
	BF16& operator= (const uint16_t& bf) noexcept {
		data = bf;
		example = *this;
		return *this;
	}
	BF16(double d) noexcept {
		example = static_cast<float>(d);
		data = reinterpret_cast<uint16_t*>(&example)[1];
	}
	operator float() const noexcept {
		uint32_t tmp = 0;
		uint16_t* tmp2 = reinterpret_cast<uint16_t*>(&tmp);
		tmp2[1] = data;
		return reinterpret_cast<float*>(&tmp)[0];
	}
	BF16(uint16_t t) noexcept {
		data = t;
		example = *this;
	}
	void print() const noexcept {
		cout << bitset<16>(data) << endl;
	}

	BF16 add(const BF16 l, const BF16 r) /*noexcept*/ {
		BF16 res;
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

	BF16 mul(const BF16 l, const BF16 r) noexcept {
		BF16 res;
		res.example = l.example * r.example;
		BF16 bfinf;
		BF16 bfnan;
		bfinf.data = 0x7f80;
		bfnan.data = 0x7f80 + 1;
		bfinf.example = INFINITY;
		bfnan.example = NAN;
		bool too_low = false;
		if (l.isbfinf() || r.isbfinf()) {
			return bfinf;
		}
		if (l.isbfnan() || r.isbfnan()) {
			return bfnan;
		}

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

	BF16 mul2(const BF16 l, const BF16 r) const noexcept {
		BF16 res;
		res.data = 0;
		res.example = l.example * r.example;
		BF16 bfinf;
		BF16 bfnan;
		bfinf.data = 0x7f80 + ((l.getsign() ^ r.getsign()) << 15);
		bfnan.data = 0x7f81 + ((l.getsign() ^ r.getsign()) << 15);
		bfinf.example = INFINITY;
		if (l.getsign() ^ r.getsign()) bfinf.example = -INFINITY;
		bfnan.example = NAN;

		if (l.isbfinf() || r.isbfinf()) {
			return bfinf;
		}
		if (l.isbfnan() || r.isbfnan()) {
			return bfnan;
		}

		uint16_t mres = ((l.getmantissa() + (uint16_t(1) << 7)) * (r.getmantissa() + (uint16_t(1) << 7)) >> 7) - (uint16_t(1) << 7);
		int16_t eres = l.getexp() + r.getexp() - int16_t(127);

		uint16_t shiftres = mres >> 7;
		
		if (eres > 0) {
			res.data += eres << 7;
			if (shiftres) {
				res.data += (mres & 0xFF80);
				res.data += (mres & 0x007F) >> shiftres;
			}
			else {
				res.data += mres;
			}
		}
		else if (eres >= -7) {
			res.data += mres;
			res.data >>= -eres; 
			res.data += uint16_t(1) << (eres + 7);
			res.data >>= 1; // +1 as subnormals
		}
		else {
			res.data = 0;
			return res;
		}
		if (res.data >> 15) {
			return bfinf;
		}
		res.data += (l.getsign() ^ r.getsign()) << 15;

		return res;
	}

	BF16& operator+= (const BF16& bf) noexcept {
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
	BF16& operator-= (const BF16& bf) noexcept {
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
	BF16& operator*= (const BF16& bf) noexcept {

	}
	BF16& operator/= (const BF16& bf) noexcept {

	}
	BF16 operator+ (const BF16& bf) const noexcept {
		BF16 res = *this;
		res += bf;
		return res;
	}
	BF16 operator- (const BF16& bf) const noexcept {
		BF16 res = *this;
		res -= bf;
		return res;
	}
	BF16 operator* (const BF16& bf) const noexcept {
		BF16 res = *this;
		res *= bf;
		return res;
	}
	BF16 operator/ (const BF16& bf) const noexcept {
		BF16 res = *this;
		res /= bf;
		return res;
	}
};

class Alltests {
	BF16 l;
	BF16 r;
public:
	Alltests() = default;
	bool equal_prec(uint16_t a, uint16_t b, uint16_t prec) {
		for (uint16_t p = 0; p <= prec; ++p) {
			if (a - p == b || a + p == b) return true;
		}
		return false;
	}
	void run() {
		uint16_t lc, rc;
		for (lc = 0x0080; lc < 0x7FFF; lc += 1) {
		//#pragma omp parallel for
			for (rc = 0x0080; rc < 0x7FFF; rc += 7) {
				cout << lc << ", " << rc;
				l = lc;
				r = rc;
				//l = l.add(l, r);
				l = l.mul2(l, r);
				//cout << endl << float(BF16(l)) << " " << float(BF16(r)) << endl;
				//cout << l.example << " " << float(l) << endl;
				if (/*!isnan(l.example)*/ (l.example == l.example) && !equal_prec(BF16(l.example).data, l.data, 1)) {
					//cout << " ADD ERROR\n";
					cout << " MUL ERROR\n";
					cout << l.example << " expected, " << float(l) << " instead\n";
					BF16(l.example).print();
					l.print();
					return;
				}
				cout << " PASSED\n";
			}
		}
	}
};

int main() {
	//uint16_t v;
	//v = 0x407f;
	//float a = BF16(v);
	//BF16(v).print();
	//v = 0xc001;
	//float b = BF16(v);
	//BF16(v).print();
	//
	//BF16 x = a + b;
	//x.print();
	//
	//cout << a << " " << b << endl;
	//cout << a + b << endl;
	// 
	//BF16 a = uint16_t(1);
	//BF16 b = uint16_t(128);
	//cout << float(a) << endl;
	//cout << float(b) << endl;
	//a.add(a,b).print();
	//cout << a.add(a,b).example << endl;

	Alltests tests;
	tests.run();

	return 0;
}