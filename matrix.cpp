#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <algorithm>

const long double PI = 3.14159265358979323846264338327950288419716939937510;


template <size_t M, size_t N, typename Field = class Rational>
class Matrix {
public:
    static const int height = M;
    static const int width = N;
    std::vector<std::vector<Field> > matrix;

    Matrix<M, N, Field> Gauss_and_det(Field& ans) const {
        Matrix<M, N, Field> new_mtr = *this;
        size_t start_i = 0;
        for (size_t j = 0; j < N; ++j) {
            size_t cur_i = start_i;
            while(cur_i < M && new_mtr.matrix[cur_i][j] == 0)
                ++cur_i;
            if (cur_i == M)
                continue;
            for (size_t i = cur_i; i > start_i; --i)
                std::swap(new_mtr.matrix[i], new_mtr.matrix[i - 1]);
            Field first_elem = new_mtr.matrix[start_i][j];
            ans *= first_elem;
            for (size_t i = 0; i < N; ++i) {
                new_mtr.matrix[start_i][i] /= first_elem;
            }
            for (size_t i = 0; i < M; ++i) {
                if (i == start_i)
                    continue;
                first_elem = new_mtr.matrix[i][j];
                for (size_t k = j; k < N; ++k)
                    new_mtr.matrix[i][k] -= new_mtr.matrix[start_i][k] * first_elem;
            }
            ++start_i;
        }
        return new_mtr;
    }


public:
    Matrix<M, N, Field>() {
        matrix.resize(M, std::vector <Field> (N, 0));
        if (M == N)
            for (size_t i = 0; i < N; ++i)
                matrix[i][i] = 1;
    }
    Matrix<M, N, Field>(int x) {
        matrix.resize(M, std::vector <Field> (N, x));
    }
    Matrix<M, N, Field>(const Matrix<M, N, Field>& construct_matrix) : matrix(construct_matrix.matrix) {}
    Matrix<M, N, Field>(std::initializer_list<std::initializer_list<int> > init_lst) {
        matrix.resize(M, std::vector<Field> (N, 0));
        int pos_i = 0, pos_j = 0;
        for (auto i : init_lst) {
            pos_j = 0;
            for (auto j : i) {
                matrix[pos_i][pos_j] =  j;
                ++pos_j;
            }
            ++pos_i;
        }
    }

    Matrix<M, N, Field>& operator+=(const Matrix<M, N, Field>& mtr) {
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                matrix[i][j] += mtr.matrix[i][j];
        return *this;
    }

    Matrix<M, N, Field> operator+(const Matrix<M, N, Field>& mtr) const {
        Matrix<M, N, Field> new_mtr(*this);
        return new_mtr += mtr;
    }

    Matrix<M, N, Field>& operator-=(const Matrix<M, N, Field>& mtr) {
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                matrix[i][j] -= mtr.matrix[i][j];
        return *this;
    }

    Matrix<M, N, Field> operator-(const Matrix<M, N, Field>& mtr) const {
        Matrix<M, N, Field> new_mtr(*this);
        return new_mtr -= mtr;
    }


    Matrix<M, N, Field> operator*(const Field& mlt) const {
        Matrix new_mtr(*this);
        return new_mtr *= mlt;
    }

    Matrix<M, N, Field>& operator*=(const Field& mlt) {
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                matrix[i][j] *= mlt;
        return *this;
    }

    Field det() const {
        static_assert(M == N, "");
        Field det_matrix = 1;
        Matrix<M, N, Field> temp_mtr = Gauss_and_det(det_matrix);
        return det_matrix;
    }

    Matrix<N, M, Field> transposed() const {
        Matrix<N, M, Field> tr_matrix;
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                tr_matrix[j][i] = matrix[i][j];
        return tr_matrix;
    }

    Field trace() const {
        static_assert(M == N, "");
        Field tr = 0;
        for (size_t i = 0; i < N; ++i)
            tr += matrix[i][i];
        return tr;
    }

    Matrix<M, N, Field> inverted() const {
        static_assert(M == N, "");
        Matrix<M, N * 2, Field> new_mtr(0);
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N * 2; ++j)
                new_mtr.matrix[i][j] = (j < N) ? matrix[i][j] : (i + N == j ? 1 : 0);
        Field temp = 0;
        new_mtr = new_mtr.Gauss_and_det(temp);
        Matrix<M, N, Field> ret_mtr = *this;
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                ret_mtr.matrix[i][j] = new_mtr.matrix[i][j + N];
        return ret_mtr;
    }

    Matrix<M, M, Field>& operator*=(const Matrix<M, M, Field>& mtr1) {
        Matrix<M, M, Field> new_mtr(0);
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < M; ++j)
                for (size_t k = 0; k < M; ++k)
                    new_mtr.matrix[i][j] += matrix[i][k] * mtr1.matrix[k][j];
        *this = new_mtr;
        return *this;
}

    void invert() {
        static_assert(M == N, "");
        *this = inverted();
        return;
    }

    Field rank() const {
        Field temp = 0;
        Matrix<M, N, Field> new_mtr = Gauss_and_det(temp);
        size_t cur = 0;
        Field ans = 0;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = cur; j < N; ++j)
                if (new_mtr.matrix[i][j] != 0) {
                    cur = j;
                    break;
                }
            if (cur < N && new_mtr.matrix[i][cur] != 0) {
                ans += 1;
                ++cur;
            }
        }
        return ans;
    }

    std::vector<Field>& operator[](int pos_i) {
        return matrix[pos_i];
    }

    std::vector<Field> operator[](int pos_i) const {
        return matrix[pos_i];
    }

    std::vector<Field> getRow(int pos_i) const {
        return matrix[pos_i];
    }

    std::vector<Field> getColumn(int pos_j) const {
        std::vector<Field> ret_column(M);
        for (size_t i = 0; i < M; ++i)
            ret_column[i] = matrix[i][pos_j];
        return ret_column;
    }

};

template <size_t N, typename Field=class Rational>
using SquareMatrix = Matrix<N, N, Field>;

template <size_t M, size_t N, size_t K, typename Field = Rational>
Matrix<M, K, Field> operator*(const Matrix<M, N, Field>& mtr1, const Matrix<N, K, Field>& mtr2) {
        Matrix<M, K, Field> new_mtr(0);
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < K; ++j)
                for (size_t k = 0; k < N; ++k)
                    new_mtr.matrix[i][j] += mtr1.matrix[i][k] * mtr2.matrix[k][j];
        return new_mtr;
}



template <size_t M, size_t N, typename Field = Rational>
Matrix<M, N, Field> operator*(const Field& num, const Matrix<M, N, Field>& mtr) {
        Matrix<M, N, Field> new_mtr(mtr);
        return new_mtr *= num;
}

template <size_t M, size_t N, typename Field = Rational>
bool operator!=(const Matrix<M, N, Field>& mtr1, const Matrix<M, N, Field>& mtr2) {
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                if (mtr1.matrix[i][j] != mtr2.matrix[i][j])
                    return true;
        return false;
}












template <int N, int K>
struct PrimeHelper {
    static const bool num = (N % K != 0) && (PrimeHelper<N, K - 1>::num);
};

template <int N>
struct PrimeHelper<N, 1> {
    static const bool num = true;
};

template <int N>
struct Prime {
    static const bool num = PrimeHelper <N, static_cast<int>(sqrt(N))>::num;
};

template <size_t N>
class Residue {
private:
    int k;

    int binpow(int x, int step) const {
        if (step == 1)
            return x;
        int p = binpow(x, step / 2);
        p = (p * p) % (int)N;
        if (step % 2)
            return (p * x) % (int)N;
        return p;
    }
public:
    Residue(int x) : k((x % (int)N + (int)N) % (int)N) {}

    Residue operator+(const Residue& res) const {
        Residue<N> new_res(k);
        return new_res += res;
    }

    Residue& operator+=(const Residue& res) {
        k = ((k + res.k) % (int)N + (int)N) % (int)N;
        return *this;
    }


    Residue operator-(const Residue& res) const {
        Residue<N> new_res(k);
        return new_res -= res;
    }

    Residue& operator-=(const Residue& res) {
        k = ((k - res.k) % (int)N + (int)N) % (int)N;
        return *this;
    }

    Residue operator*(const Residue& res) const {
        Residue<N> new_res(k);
        return new_res *= res;
    }

    Residue& operator*=(const Residue& res) {
        k = ((k * res.k) % (int)N + (int)N) % (int)N;
        return *this;
    }

    Residue operator/(const Residue& res) const {
        Residue<N> new_res(k);
        return new_res /= res;
    }

    Residue operator/=(const Residue& res) {
        static_assert(Prime<N>::num, "");
        k = ((k * binpow((res.k % (int)N + (int)N) % (int)N, (int)N - 2)) % (int)N + (int)N) % (int)N;
        return *this;
    }


    explicit operator int() const {
        return k;
    }

};


template <size_t N>
bool operator==(const Residue<N>& res1, const Residue<N>& res2) {
        return res1.k == res2.k;
}

template <size_t N>
bool operator!=(const Residue<N>& res1, const Residue<N>& res2) {
        return res1.k != res2.k;
}








using namespace std;

class BigInteger{
public:
    vector <long long> num;
    int osn = 3;
    int step = 1000;
    bool sign_minus = false;



    bool operator_greater_for_abs(const BigInteger& a, const BigInteger& b) const{
        if (a.num.size() != b.num.size())
            return (a.num.size() > b.num.size());
        bool fl = false;
        bool eq = true;
        for (int i = a.num.size() - 1; i >= 0; --i){
            if (a.num[i] != b.num[i]){
                eq = false;
                if (a.num[i] > b.num[i])
                    fl = true;
                break;
            }

        }
        return (eq || fl);
    }

    bool equivalent(const BigInteger& a, const BigInteger& b) const{
        if (a.num.size() != b.num.size())
            return false;
        for (size_t i = 0; i < a.num.size(); ++i)
            if (a.num[i] != b.num[i])
                return false;
        return true;
    }

    int inv(int x, int t){
        int k = 0;
        for (int i = 0; i < t; ++i)
            if ((x & (1 << i)) != 0)
                k += (1 << (t - i - 1));
        return k;
    }


    void fft(std::vector <complex <long double> > &cur, int rev, int t, int mx){
        for (int i = 0; i < mx; ++i){
            int rv = inv(i, t);
            if (i >= rv)
                swap(cur[i], cur[rv]);
        }
        for (int len = 2; len <= mx; len *= 2){
            for (int i = 0; i < mx; i += len){
                std::complex <long double> ang(cos(rev * 2 * PI / len), sin(rev * 2 * PI / len));
                std::complex <long double> cur_ang(1, 0);
                for (int j = i; j < i + len / 2; ++j){
                    std::complex <long double> ff = cur[j];
                    std::complex <long double> ss = cur[j + len / 2];
                    cur[j] = ff + cur_ang * ss;
                    cur[j + len / 2] = ff - cur_ang * ss;
                    cur_ang *= ang;
                }
            }
        }
    }

public:

    BigInteger(){}

    BigInteger(int number){
        if (number < 0){
            sign_minus = true;
            number = -number;
        }
        if (number == 0)
            num.push_back(0);
        while(number != 0){
            num.push_back(number % 10);
            number /= 10;
        }
    }

    void construct_int(int number){
        if (number < 0){
            sign_minus = true;
            number = -number;
        }
        if (number == 0)
            num.push_back(0);
        while(number != 0){
            num.push_back(number % 10);
            number /= 10;
        }
    }

    BigInteger(const BigInteger& a){
        num = a.num;
        sign_minus = a.sign_minus;
    }

    void construct_bg(const BigInteger& a){
        num = a.num;
        sign_minus = a.sign_minus;
    }


    BigInteger& operator+=(const BigInteger& bg){
        if (((int)bg.sign_minus + (int)sign_minus) % 2){
            sign_minus ^= 1;
            *this -= bg;
            sign_minus ^= 1;
            return *this;
        }
        int fl = 0;
        size_t sz = max(num.size(), bg.num.size());
        for (size_t i = 0; i < sz; ++i){
            int k = fl;
            if (i < num.size())
                k += num[i];
            if (i < bg.num.size())
                k += bg.num[i];
            if (i < num.size())
                num[i] = k % 10;
            else
                num.push_back(k % 10);
            fl = k / 10;
        }
        if (fl)
            num.push_back(fl);
        prov_null();
        return *this;
    }



    BigInteger& operator-=(const BigInteger& bg){
        if (((int)bg.sign_minus + (int)sign_minus) % 2){
            sign_minus ^= 1;
            *this += bg;
            sign_minus ^= 1;
            return *this;
        }
        bool fl = operator_greater_for_abs(*this, bg);
        bool tr = false;
        if (!fl){
            sign_minus = (1 ^ bg.sign_minus);
            while(num.size() < bg.num.size())
                num.push_back(0);
            for (size_t i = 0; i < num.size(); ++i){
                int k = bg.num[i] - tr - num[i];
                tr = 0;
                if (k < 0)
                    k += 10, tr = 1;
                num[i] = k;
            }
        } else{
            for (size_t i = 0; i < num.size(); ++i){
                int k = num[i] - tr - (i < bg.num.size() ? bg.num[i] : 0);
                tr = 0;
                if (k < 0)
                    k += 10, tr = 1;
                num[i] = k;
            }
        }
        while(num.size() > 1 && num.back() == 0)
            num.pop_back();
        prov_null();
        return *this;
    }




    BigInteger& operator%=(const BigInteger& a){
        return (*this) -= (*this) / a * a;
    }



    BigInteger operator/(const BigInteger& b) const{
        BigInteger a(*this);
        return a /= b;
    }

    BigInteger& operator/=(const BigInteger& a){
        BigInteger pp(*this), pp1(a);
        pp1.sign_minus = pp.sign_minus = false;
        std::reverse(pp.num.begin(), pp.num.end());
        BigInteger ans, cur;
        ans.num.push_back(0);
        bool fl = false;
        BigInteger x(pp1);
        for (size_t i = 0; i < pp.num.size(); ++i){
            if (cur.num.size() == 1 && cur.num[0] == 0)
                cur.num.pop_back();
            cur.num.push_back(pp.num[i]);
            std::reverse(cur.num.begin(), cur.num.end());
            if (operator_greater_for_abs(cur, pp1)){
                fl = true;
                for (int j = 2; j <= 10; ++j){
                    x += pp1;
                    if (equivalent(x, cur)){
                        ans.num.push_back(j);
                        cur -= x;
                        for (int k = 0; k < j - 1; ++k)
                            x -= pp1;
                        break;
                    }
                    if (j == 10 || operator_greater_for_abs(x, cur)){
                        ans.num.push_back(j - 1);
                        cur -= (x -= pp1);
                        for (int k = 0; k < j - 2; ++k)
                            x -= pp1;
                        break;
                    }
                }
            }else{
                if (fl)
                    ans.num.push_back(0);
            }
            std::reverse(cur.num.begin(), cur.num.end());
        }
        std::reverse(ans.num.begin(), ans.num.end());
        while(ans.num.size() > 1 && ans.num.back() == 0)
            ans.num.pop_back();
        num.resize(ans.num.size());
        for (size_t i = 0; i < ans.num.size(); ++i){
            num[i] = ans.num[i];
        }
        sign_minus ^= a.sign_minus;
        prov_null();
        return *this;
    }

    std::string division(const BigInteger& a, int k) const{
        BigInteger pp(*this), pp1(a);
        std::string s = "";
        if ((*this) && (sign_minus ^ a.sign_minus))
            s += '-';
        pp1.sign_minus = pp.sign_minus = false;
        std::reverse(pp.num.begin(), pp.num.end());
        BigInteger ans, cur;
        ans.num.push_back(0);
        bool fl = false;
        int t = 0;
        BigInteger x(pp1);
        for (size_t i = 0; i < pp.num.size() || t <= k; ++i){
            if (cur.num.size() == 1 && cur.num[0] == 0)
                cur.num.pop_back();
            if (i < pp.num.size())
                cur.num.push_back(pp.num[i]);
            else{
                ++t;
                cur.num.push_back(0);
            }
            std::reverse(cur.num.begin(), cur.num.end());
            if (operator_greater_for_abs(cur, pp1)){
                fl = true;
                for (int j = 2; j <= 10; ++j){
                    x += pp1;
                    if (equivalent(x, cur)){
                        cur -= x;
                        for (int k = 0; k < j - 1; ++k)
                            x -= pp1;
                        s += char(j + '0');
                        break;
                    }
                    if (j == 10 || operator_greater_for_abs(x, cur)){
                        s += char(j - 1 + '0');
                        cur -= (x -= pp1);
                        for (int k = 0; k < j - 2; ++k)
                                x -= pp1;
                        break;
                    }
                }
            }else{
                s += '0';
            }
            std::reverse(cur.num.begin(), cur.num.end());
        }
        fl = 0;
        if (s.back() - '0' >= 5)
            fl = 1;
        s.pop_back();
        std::string s1 = "";
        while(s.size()){
            if ((int)s1.size() == k && k != 0)
                s1 += '.';
            if (s.back() == '-'){
                s1 += s.back();
                s.pop_back();
                continue;
            }
            int k = s.back() - '0' + fl;
            s1 += char(k % 10 + '0');
            fl = k / 10;
            s.pop_back();
        }
        if (s1.back() == '-'){
            s += '-';
            s1.pop_back();
        }
        while(s1.size() > 1 && s1.back() != '.' && s1.back() == '0')
            s1.pop_back();
        if (s1.size() > 0 && s1.back() == '.')
            s1 += '0';
        for (int i = s1.length() - 1; i >= 0; --i){
            s += s1[i];
        }
        return s;
    }



    BigInteger operator*(const BigInteger& b) const{
        BigInteger a(*this);
        return a *= b;
    }

    BigInteger& operator*=(const BigInteger& a){
        std::vector <complex <long double> > pp, pp1;
        int k = 0;
        std::complex <long double> x(0, 0);
        for (size_t i = 0; i < a.num.size(); ++i){
            ++k;
            std::complex <long double> t(a.num[i], 0);
            int tt = 1;
            for (size_t j = 0; j < (size_t)(((k % osn) - 1 + osn) % osn); ++j)
                tt *= 10;
            std::complex <long double> xx(tt, 0);
            x += xx * t;
            if (k % osn == 0 || i == a.num.size() - 1){

                pp.push_back(x);
                std::complex <long double> xz(0, 0);
                x = xz;
            }
        }
        k = 0;
        std::complex <long double> xz(0, 0);
        x = xz;
        for (size_t i = 0; i < num.size(); ++i){
            ++k;
            std::complex <long double> t(num[i], 0);
            int tt = 1;
            for (size_t j = 0; j < (size_t)(((k % osn) - 1 + osn) % osn); ++j)
                tt *= 10;
            std::complex <long double> xx(tt, 0);
            x += xx * t;
            if (k % osn == 0 || i == num.size() - 1){
                pp1.push_back(x);
                std::complex <long double> xz(0, 0);
                x = xz;
            }
        }
        int mx = max(pp.size(), pp1.size());
        mx *= 2;
        int t = 0;
        while((1 << t) < mx)
            ++t;
        mx = (1 << t);
        while(pp.size() != (size_t)mx){
            std::complex <long double> x(0, 0);
            pp.push_back(x);
        }
        while(pp1.size() != (size_t)mx){
            std::complex <long double> x(0, 0);
            pp1.push_back(x);
        }
        fft(pp, 1, t, mx);
        fft(pp1, 1, t, mx);
        for (int i = 0; i < mx; ++i)
            pp[i] *= pp1[i];
        fft(pp, -1, t, mx);
        std::complex <long double> xx(mx, 0);
        x = xx;
        std::vector <char> ans;
        long long fl = 0;
        for (size_t i = 0; i < (size_t)mx; ++i){
            long long t = round((pp[i] / x).real());
            fl += t;
            string s = "";
            int k = fl % step;
            while(k != 0){
                s += char('0' + k % 10);
                k /= 10;
            }
            while(s.length() != (size_t)osn)
                s += '0';
            for (size_t j = 0; j < s.length(); ++j)
                ans.push_back(s[j]);
            fl /= step;
        }
        while(ans.size() > 1 && ans.back() == '0')
            ans.pop_back();
        num.resize(ans.size());
        for (size_t i = 0; i < ans.size(); ++i){
            num[i] = ans[i] - '0';
        }
        sign_minus ^= a.sign_minus;
        prov_null();
        return *this;
    }




    bool operator<(const BigInteger& a) const{
        if (sign_minus && !a.sign_minus)
            return true;
        if (!sign_minus && a.sign_minus)
            return false;
        if (num.size() != a.num.size())
            return (num.size() < a.num.size()) ^ sign_minus;
        bool fl = true;
        bool eq = true;
        for (int i = num.size() - 1; i >= 0; --i){
            if (num[i] != a.num[i]){
                eq = false;
                if (num[i] > a.num[i])
                    fl = false;
                break;
            }

        }
        if (eq)
            return false;
        return fl ^ sign_minus;
    }


    BigInteger& operator++(){
        if (sign_minus){
            sign_minus = false;
            --(*this);
            sign_minus = true;
            prov_null();
            return *this;
        }
        int fl = 1;
        for (size_t i = 0; i < num.size(); ++i){
            int k = num[i] + fl;
            num[i] = k % 10;
            fl = k / 10;
            if (!fl)
                break;
        }
        if (fl)
            num.push_back(1);
        prov_null();
        return *this;
    }

    BigInteger& operator--(){
        if (sign_minus){
            sign_minus = false;
            ++(*this);
            sign_minus = true;
            prov_null();
            return *this;
        }
        if (!(*this)){
            num[0] = 1;
            sign_minus = true;
            return *this;
        }
        for (size_t i = 0; i < num.size(); ++i){
            if (num[i] != 0){
                num[i] -= 1;
                break;
            }
            num[i] = 9;
        }
        if (num.size() > 1 && num.back() == 0)
            num.pop_back();
        prov_null();
        return *this;
    }

    BigInteger operator-() const{
        BigInteger a(*this);
        a.sign_minus ^= true;
        a.prov_null();
        return a;
    }

    BigInteger operator++(int) {
        if (sign_minus){
            BigInteger x(*this);
            sign_minus = false;
            --(*this);
            sign_minus = true;
            prov_null();
            return x;
        }
        BigInteger x(*this);
        ++(*this);
        return x;
    }

    BigInteger operator--(int) {
        if (sign_minus){
            BigInteger x(*this);
            sign_minus = false;
            ++(*this);
            sign_minus = true;
            prov_null();
            return x;
        }
        BigInteger x(*this);
        ++(*this);
        return x;
    }

    explicit operator bool() const{
        return !(num.size() == 1 && num[0] == 0);
    }

    explicit operator int() const{
        int a = 0;
        for (int i = num.size() - 1; i >= 0; --i)
            a = a * 10 + num[i];
        if (sign_minus)
            a *= -1;
        return a;
    }

    std::string toString() const{
        std::string s = "";
        if (sign_minus)
            s += '-';
        for (int i = num.size() - 1; i >= 0; --i){
            s += char(num[i] + '0');
        }
        return s;
    }


    friend std::ostream& operator<<(std::ostream& out, const BigInteger& a){
        if (a.sign_minus)
            out << '-';
        for (int i = a.num.size() - 1; i >= 0; --i)
            out << a.num[i];
        return out;
    }

    friend std::istream& operator>>(std::istream& in, BigInteger& a){
        a.clear();
        char c;
        while(true){
            c = in.get();
            if (!isspace(c)){
                if (c == '-')
                    a.sign_minus = true;
                else
                    a.num.push_back(c - '0');
                break;
            }
        }
        while (true){
            c = in.get();
            if (isspace(c) || c == EOF)
                break;
            a.num.push_back(c - '0');
        }
        std::reverse(a.num.begin(), a.num.end());
        return in;
    }

    void clear(){
        num.resize(0);
        sign_minus = false;
    }

    void prov_null() {
        if (num.size() == 1 && num[0] == 0)
            sign_minus = false;
    }



    void change_sign_minus(bool x) {
        sign_minus ^= x;
    }

    void assign_sign_minus(bool x) {
        sign_minus = x;
    }

    bool get_sign_minus() const {
        return sign_minus;
    }

    bool check_null() const {
        return (num.size() == 1 && num[0] == 0);
    }

    bool IsEven() const {
        return !(num[0] % 2);
    }

    void divide_two() {
        int BaseOverflow = 0;
        for (int i = num.size() - 1; i >= 0; --i) {
            int k = num[i] + BaseOverflow * 10;
            num[i] = k / 2;
            BaseOverflow = k % 2;
        }
        while(num.size() > 1 && num.back() == 0)
            num.pop_back();
    }

    void multiply_two() {
        int BaseOverflow = 0;
        for (size_t i = 0; i < num.size(); ++i) {
            num[i] = num[i] * 2 + BaseOverflow;
            BaseOverflow = num[i] / 10;
            num[i] %= 10;
        }
        if (BaseOverflow)
            num.push_back(1);
    }

    ~BigInteger(){
        num.resize(0);
    };

};

BigInteger operator-(const BigInteger& a, const BigInteger& b){
    BigInteger x(a);
    return x -= b;
}

BigInteger operator+(const BigInteger& a, const BigInteger& b){
    BigInteger x(a);
    return x += b;
}


BigInteger operator*(int a, const BigInteger& b){
    BigInteger x(a);
    return x *= b;
}

BigInteger operator/(int a, const BigInteger& b){
    BigInteger x(a);
    return x /= b;
}

BigInteger operator%(const BigInteger& a, const BigInteger& b){
    BigInteger x(a);
    return x %= b;
}

bool operator<(int a, const BigInteger& b){
    BigInteger a1(a);
    return a1 < b;
}

bool operator>(const BigInteger& a, const BigInteger& b){
    return -a < -b;
}

bool operator<=(const BigInteger& a, const BigInteger& b){
    return !(b < a);
}

bool operator>=(const BigInteger& a, const BigInteger& b){
    return !(a < b);
}

bool operator!=(const BigInteger& a, const BigInteger& b){
    return (a < b || b < a);
}

bool operator==(const BigInteger& a, const BigInteger& b){
    return !(a < b || b < a);
}


class Rational{
private:
    BigInteger numerator, denominator;
public:

    Rational(){}

    Rational(int a){
        numerator = a;
        denominator = 1;
    }

    Rational(const BigInteger& a){
        numerator = a;
        denominator = 1;
    }

    Rational(const Rational& a): numerator(a.numerator), denominator(a.denominator){}

    BigInteger gcd(BigInteger a, BigInteger b) const{
        BigInteger ans = 1;
        while (!a.check_null() && !b.check_null()) {
            int cnt_even = 0;
            if (a.IsEven()) {
                a.divide_two();
                ++cnt_even;
            }
            if (b.IsEven()) {
                b.divide_two();
                ++cnt_even;
            }
            if (cnt_even == 2)
                ans.multiply_two();
            if (!a.IsEven() && !b.IsEven()) {
                if (a > b)
                    a -= b;
                else
                    b -= a;
            }
        }

        return ans *= ((b == 0) ? a : b);
    }

    void reduce(){
        numerator.change_sign_minus(denominator.get_sign_minus());
        denominator.change_sign_minus(denominator.get_sign_minus());
        bool fl = numerator.get_sign_minus();
        numerator.assign_sign_minus(false);
        BigInteger x(gcd(denominator, numerator));
        numerator.change_sign_minus(fl);
        numerator /= x;
        denominator /= x;
    }



    Rational& operator+=(const Rational& a){
        numerator.change_sign_minus(denominator.get_sign_minus());
        denominator.change_sign_minus(denominator.get_sign_minus());
        numerator *= a.denominator;
        numerator += a.numerator * denominator;
        denominator *= a.denominator;
        reduce();
        return *this;
    }



    Rational& operator-=(const Rational& a){
        numerator.change_sign_minus(denominator.get_sign_minus());
        denominator.change_sign_minus(denominator.get_sign_minus());
        numerator *= a.denominator;
        numerator -= a.numerator * denominator;
        denominator *= a.denominator;
        reduce();
        return *this;
    }


    Rational& operator*=(const Rational& a){
        numerator.change_sign_minus(denominator.get_sign_minus());
        denominator.change_sign_minus(denominator.get_sign_minus());
        numerator *= a.numerator;
        denominator *= a.denominator;
        reduce();
        return *this;
    }



    Rational& operator/=(const Rational& a){
        numerator.change_sign_minus(denominator.get_sign_minus());
        denominator.change_sign_minus(denominator.get_sign_minus());
        numerator *= a.denominator;
        denominator *= a.numerator;
        reduce();
        return *this;
    }


    Rational operator-() const{
        Rational a(*this);
        a.numerator.change_sign_minus(true);
        return a;
    }

    explicit operator double() const{
        string s = asDecimal(24);
        long double ans = 0;
        long double k = 1;
        size_t start = (size_t)(s[0] == '-');
        bool fl = 0;
        for (size_t i = start; i < s.length(); ++i){
            if (s[i] == '.'){
                fl = true;
                continue;
            }
            if (!fl)
                ans = ans * 10 + s[i] - '0';
            else{
                k *= 10;
                ans += (s[i] - '0') * (1.0 / k);
            }
        }
        if (s[0] == '-')
            ans *= -1;
        return ans;
    }

    explicit operator int() const{
        int a = (int)numerator;
        if (numerator.sign_minus)
            a *= -1;
        return a;
    }



    bool operator<(const Rational& a) const{
        if ((numerator.get_sign_minus() ^ denominator.get_sign_minus()) && !(a.numerator.get_sign_minus() ^ a.denominator.get_sign_minus()))
            return true;
        if (!(numerator.get_sign_minus() ^ denominator.get_sign_minus()) && (a.numerator.get_sign_minus() ^ a.denominator.get_sign_minus()))
            return false;
        BigInteger x;
        x = numerator * a.denominator;
        x.assign_sign_minus(numerator.get_sign_minus() ^ denominator.get_sign_minus());
        if (a.denominator.get_sign_minus())
            x.change_sign_minus(true);
        BigInteger x1;
        x1 = a.numerator * denominator;
        x1.assign_sign_minus(a.numerator.get_sign_minus() ^ a.denominator.get_sign_minus());
        if (denominator.get_sign_minus())
            x1.change_sign_minus(true);
        return x < x1;
    }

    string toString() {
        reduce();
        string s = "";
        numerator.assign_sign_minus(numerator.get_sign_minus() ^ denominator.get_sign_minus());
        denominator.assign_sign_minus(false);
        s += numerator.toString();
        if (!numerator || !(denominator - 1))
            return s;
        s += '/';
        s += denominator.toString();
        return s;
    }

    string asDecimal(int k = 0) const{
        return numerator.division(denominator, k);
    }


    friend std::istream& operator>>(std::istream& in, Rational& a){
        a.denominator.clear();
        a.denominator = 1;
        a.numerator.clear();
        int new_num = 0;
        char c;
        bool isneg = false;
        while(true){
            c = in.get();
            if (!isspace(c)){
                if (c == '-')
                    isneg = 1;
                else
                    new_num = new_num * 10 + (c - '0');
                break;
            }
        }
        while (true){
            c = in.get();
            if (isspace(c) || c == EOF)
                break;
            new_num = new_num * 10 + (c - '0');
        }
        if (isneg)
            new_num *= -1;
        a.numerator = new_num;
        //cout << a.numerator.sign_minus << endl;
        return in;
    }

};



Rational operator-(const Rational& a, const Rational& b){
    Rational x(a);
    return x -= b;
}

Rational operator+(const Rational& a, const Rational& b){
    Rational x(a);
    return x += b;
}


Rational operator*(const Rational& a, const Rational& b){
    Rational x(a);
    return x *= b;
}


Rational operator/(const Rational& a, const Rational& b){
    Rational x(a);
    return x /= b;
}


bool operator<(int a, const Rational& b){
    Rational a1(a);
    return a1 < b;
}

bool operator>(const Rational& a, const Rational& b){
    return -a < -b;
}

bool operator<=(const Rational& a, const Rational& b){
    return !(b < a);
}

bool operator>=(const Rational& a, const Rational& b){
    return !(a < b);
}

bool operator!=(const Rational& a, const Rational& b){
    return (a < b || b < a);
}

bool operator==(const Rational& a, const Rational& b){
    return !(a < b || b < a);
}

