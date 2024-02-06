
#pragma once

#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>



using namespace std;

// Utilisation de la classe mean_var définie au cours 

// Moyenne, variance, IC (MC)
struct mean_var {
    mean_var(double n = 0, double sum_x = 0, double sum_xx = 0)
        : sample_size(n), sum_x(sum_x), sum_xx(sum_xx) { }
    [[nodiscard]] double mean() const { return sum_x / (double) sample_size; }
    [[nodiscard]] double var() const { return (sum_xx - sample_size * mean() * mean())
        / (double) (sample_size-1); }
    [[nodiscard]] double ic_size() const { return 1.96 * std::sqrt(var() / sample_size); }

    mean_var & operator+=(mean_var const & mv);
    friend mean_var operator+(mean_var const & mv1, mean_var const & mv2);
    friend std::ostream & operator<<(std::ostream & o, mean_var const & mv);

    private:
        double sample_size;
        double sum_x, sum_xx;
};

// definition of the method operator+=
mean_var & mean_var::operator+=(mean_var const & mv) {
    sample_size += mv.sample_size;
    sum_x += mv.sum_x;
    sum_xx += mv.sum_xx;
    return *this;
}

// definition of the global function operator+ declared friend in the class mean_var
mean_var operator+(mean_var const & mv1, mean_var const & mv2) { //permet de paralléliser!
    return { mv1.sample_size + mv2.sample_size,
        mv1.sum_x + mv2.sum_x,
        mv1.sum_xx + mv2.sum_xx };
}

std::ostream & operator<<(std::ostream & o, mean_var const & mv) {
    return o << "Size: " << mv.sample_size
        << "\tMean: " << mv.mean()
        << "\tVar: " << mv.var();
}


vector<double> normal_generation(unsigned size){
    std::normal_distribution<double> X; 
    random_device rd;
    auto seed = rd();
    mt19937_64 gen(seed);
    vector<double> brow(size);
    generate(brow.begin(), brow.end(), [&](){return X(gen);});
    return brow;
}
vector<double> negative_tab(const vector<double>& X){
    vector<double> out(X.size());
    transform(cbegin(X), cend(X), begin(out), [] (double x) {return -x;});
    return out;
}