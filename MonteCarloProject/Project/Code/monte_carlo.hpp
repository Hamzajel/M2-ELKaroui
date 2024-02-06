#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <fstream>

#include "mc.hpp"
#include "payoff.hpp"
#include "diff_finie.hpp"
#include "tangent_process.hpp"
#include "malliavin.hpp"
#include "pseudo_CEV.hpp"
#include "BS.hpp"

using namespace std;

//Fonction qui permet de générer les calculs d'espérance Monte Carlo

template <typename T_method>
mean_var monte_carlo(T_method* obj, double (T_method::*fp)(string, vector<double>), double sample_size, string payoff, unsigned size_h) {
    vector<double> brow;
    double sum_x = 0, sum_xx = 0;
    for (unsigned k = 0; k < sample_size; ++k) {
        brow=normal_generation(size_h);
        double x = (obj->*fp)(payoff, brow);
        sum_x += x;
        sum_xx += x*x;
    }
    return { sample_size, sum_x, sum_xx };
}

template <typename T_method>
mean_var monte_carlo_anti(T_method* obj, double (T_method::*fp)(string, vector<double>), double sample_size, string payoff, unsigned size_h) {
    vector<double> brow;
    double sum_x = 0, sum_xx = 0;
    for (unsigned k = 0; k < sample_size; ++k) {
        brow=normal_generation(size_h);
        double x = 0.5*((obj->*fp)(payoff, brow)+(obj->*fp)(payoff, negative_tab(brow)));
        sum_x += x;
        sum_xx += x*x;
    }
    return { sample_size, sum_x, sum_xx };
}

template <typename T_method> // attention que pour le delta et Call/Put
mean_var monte_carlo_con_var(T_method* obj, double (T_method::*fp)(string, vector<double>, double), unsigned sample_size, string payoff, unsigned size_h) {
    vector<double> brow;
    double sum_x = 0, sum_xx = 0;
    auto size_L=0.1*sample_size;
    double lambda=obj->lambda_delta(payoff,size_L);

    auto size_MC=0.9*sample_size;
    for (unsigned k = 0; k < size_MC; ++k) {
        brow=normal_generation(size_h);
        double x = (obj->*fp)(payoff, brow, lambda);
        sum_x += x;
        sum_xx += x*x;
    }
    return { size_MC, sum_x, sum_xx };
}

template <typename T_method> // attention que pour le delta et Call/Put
mean_var monte_carlo_con_var_anti(T_method* obj, double (T_method::*fp)(string, vector<double>, double), unsigned sample_size, string payoff, unsigned size_h) {
    vector<double> brow;
    double sum_x = 0, sum_xx = 0;
    auto size_L=0.1*sample_size;
    double lambda=obj->lambda_delta(payoff,size_L);

    auto size_MC=0.9*sample_size;
    for (unsigned k = 0; k < size_MC; ++k) {
        brow=normal_generation(size_h);
        double x = 0.5*((obj->*fp)(payoff, brow, lambda)+(obj->*fp)(payoff, negative_tab(brow), lambda));
        sum_x += x;
        sum_xx += x*x;
    }
    return { size_MC, sum_x, sum_xx };
}