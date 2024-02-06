#include <iostream>
#include <random>
#include "monte_carlo.hpp"

using namespace std;




int main(){

    random_device rd;
    auto seed = rd();
    mt19937_64 gen(seed);
    double S=100,r=0.1, sigma=0.6,eps=1e-2,T=1,K=100;
    string type_payoff="call";
    auto sample_size_max=1e5; // M=o(eps^-4) si payoff regular !!
    auto sample_size_init = 10000;
    auto step = 1000;
    auto n_cev=1;
    double v=0.9;

    std::cout<<"Sample size B&S: "<<sample_size_max<<std::endl;

    std::cout<<"-------------------------------            "<< type_payoff<<"             -----------------------------------"<<std::endl;

    std::cout<<"-------------------------------                   B&S                     -----------------------------------"<<std::endl;

    // Diff finie à pas constant 
   simu_payoff_BS sim(S,r,sigma,T,K);

   

    diff_finie<simu_payoff_BS> df(sim,eps,1,r,T);

    std::cout<<"-------------------------------    Difference finie, pas constant, B&S   -----------------------------------"<<std::endl;

    std::ofstream outfile;
    outfile.open("Diff_finie_cst_BS.csv");
    outfile << "Sample Size, Greek, IC, Greek Anti, IC Anti, Greek Con Var, IC Con Var\n";

    for (unsigned sample_size = sample_size_init; sample_size <= sample_size_max; sample_size += step) {
        auto greek = monte_carlo(&df, &diff_finie<simu_payoff_BS>::delta, sample_size, type_payoff, 1);
        auto greek_anti = monte_carlo_anti(&df, &diff_finie<simu_payoff_BS>::delta, sample_size, type_payoff, 1);
        auto greek_var_con = monte_carlo_con_var(&df, &diff_finie<simu_payoff_BS>::delta_con_var, sample_size, type_payoff, 1);
        auto greek_var_con_anti = monte_carlo_con_var_anti(&df, &diff_finie<simu_payoff_BS>::delta_con_var, sample_size, type_payoff, 1);

        outfile << sample_size << ",";
        outfile << greek.mean() << ",";
        outfile << greek.ic_size() << ",";
        outfile << greek_anti.mean() << ",";
        outfile << greek_anti.ic_size() << ",";
        outfile << greek_var_con.mean() << ",";
        outfile << greek_var_con.ic_size() << ",";
        outfile << greek_var_con_anti.mean() << ",";
        outfile << greek_var_con_anti.ic_size() << "\n";

    }

    outfile.close();
    // Output a message to the console indicating that the CSV file was written
    std::cout << "Results written to Diff_finie_cst_BS.csv" << std::endl;

    /*cout<<"greek "<< greek.mean()<<" IC: "<<greek.ic_size()<<endl;
    cout<<"greek anti "<< greek_anti.mean()<<" IC: "<<greek_anti.ic_size()<<endl;
    cout<<"greek con_var "<< greek_var_con.mean()<<" IC: "<<greek_var_con.ic_size()<<endl;
    cout<<"greek con_var_anti "<< greek_var_con_anti.mean()<<" IC: "<<greek_var_con_anti.ic_size()<<endl;*/

    std::cout<<"-------------------------------    Difference finie, pas decroissant, B&S   -----------------------------------"<<std::endl;
    double init_eps=1e-3;

    std::ofstream outfile2;
    outfile2.open("Diff_finie_dec_BS.csv");
    outfile2 << "Sample Size, Greek, IC, Greek Anti, IC Anti, Greek Con Var, IC Con Var\n";

    for (unsigned sample_size = sample_size_init; sample_size <= sample_size_max; sample_size += step) {
        auto greek_df_d = df.MC_decrease(sample_size, type_payoff, init_eps, 1);
        auto greek_df_d_anti = df.MC_decrease_anti(sample_size, type_payoff, init_eps, 1);
        auto greek_df_d_con_var = df.MC_decrease_cont_var(r, T, sample_size, type_payoff, init_eps, 1);
        auto greek_var_con_anti = monte_carlo_con_var_anti(&df,&diff_finie<simu_payoff_BS>::delta_con_var,sample_size,type_payoff, 1);

        outfile2 << greek_df_d[2].mean() << ",";
        outfile2 << greek_df_d[0].ic_size() << ",";
        outfile2 << greek_df_d_anti[2].mean() << ",";
        outfile2 << greek_df_d_anti[0].ic_size() << ",";
        outfile2 << greek_df_d_con_var.mean() << ",";
        outfile2 << greek_df_d_con_var.ic_size() << ",";
        outfile2 << greek_var_con_anti.mean() << ",";
        outfile2<< greek_var_con_anti.ic_size() << std::endl;

        /*cout << "greek " << greek_df_d[2].mean() << " IC: " << greek_df_d[0].ic_size() << endl;
        cout << "greek anti " << greek_df_d_anti[2].mean() << " IC: " << greek_df_d_anti[0].ic_size() << endl;
        cout << "greek con_var " << greek_df_d_con_var.mean() << " IC: " << greek_df_d_con_var.ic_size() << endl;
        cout << "greek con_var_anti " << greek_var_con_anti.mean() << " IC: " << greek_var_con_anti.ic_size() << endl;*/
    }

    outfile2.close();
    // Output a message to the console indicating that the CSV file was written
    std::cout << "Results written to Diff_finie_dec_BS.csv" << std::endl;


   // Flot tangent 
    tangent_process<simu_payoff_BS> tan(r,sample_size_max,T,K,sim,1);

    std::ofstream outfile3;
    outfile3.open("Flot_tangent_BS.csv");
    outfile3 << "Sample Size, Greek, IC, Greek Anti, IC Anti, Greek Con Var, IC Con Var, Greek Anti Con Var, IC Anti Con Var\n";


    std::cout<<"-------------------------------            Flot tangent, B&S             -----------------------------------"<<std::endl;
    for (unsigned sample_size = sample_size_init; sample_size <= sample_size_max ; sample_size+=step) {
        auto greek = monte_carlo(&tan, &tangent_process<simu_payoff_BS>::delta, sample_size, type_payoff, 1);
        auto greek_anti = monte_carlo_anti(&tan, &tangent_process<simu_payoff_BS>::delta, sample_size, type_payoff, 1);
        auto greek_var_con = monte_carlo_con_var(&tan, &tangent_process<simu_payoff_BS>::delta_con_var, sample_size,
                                            type_payoff, 1);
        auto greek_var_con_anti = monte_carlo_con_var_anti(&tan, &tangent_process<simu_payoff_BS>::delta_con_var,
                                                      sample_size, type_payoff, 1);

        // Write the calculated values to the CSV file
        outfile3 << sample_size << ",";
        outfile3 << greek.mean() << ",";
        outfile3 << greek.ic_size() << ",";
        outfile3 << greek_anti.mean() << ",";
        outfile3 << greek_anti.ic_size() << ",";
        outfile3 << greek_var_con.mean() << ",";
        outfile3 << greek_var_con.ic_size() << ",";
        outfile3 << greek_var_con_anti.mean() << ",";
        outfile3 << greek_var_con_anti.ic_size() << "\n";

        /*cout << "greek " << greek.mean() << " IC: " << greek.ic_size() << endl;
        cout << "greek anti " << greek_anti.mean() << " IC: " << greek_anti.ic_size() << endl;
        cout << "greek con_var " << greek_var_con.mean() << " IC: " << greek_var_con.ic_size() << endl;
        cout << "greek con_var_anti " << greek_var_con_anti.mean() << " IC: " << greek_var_con_anti.ic_size() << endl;*/
    }

    // Close the CSV file
    outfile3.close();
    // Output a message to the console indicating that the CSV file was written
    std::cout << "Results written to Flot_tangent_BS.csv" << std::endl;

    malliavin<simu_payoff_BS> mal(r,sample_size_max,T,K,sim,1);

    std::cout<<"-------------------------------             Malliavin, B&S              -----------------------------------"<<std::endl;

    std::ofstream outfile4;
    outfile4.open("Malliavin_BS.csv");
    outfile4 << "Sample Size, Greek, IC, Greek Anti, IC Anti, Greek Con Var, IC Con Var\n";

    for (unsigned sample_size = sample_size_init; sample_size <= sample_size_max ; sample_size+=step) {
        auto greek = monte_carlo(&mal, &malliavin<simu_payoff_BS>::delta, sample_size, type_payoff, 1);
        auto greek_anti = monte_carlo_anti(&mal, &malliavin<simu_payoff_BS>::delta, sample_size, type_payoff, 1);
        auto greek_var_con = monte_carlo_con_var(&mal, &malliavin<simu_payoff_BS>::delta_con_var, sample_size, type_payoff,
                                            1);
        auto greek_var_con_anti = monte_carlo_con_var_anti(&mal, &malliavin<simu_payoff_BS>::delta_con_var, sample_size,
                                                      type_payoff, 1);

        outfile4 << sample_size << ",";
        outfile4 << greek.mean() << ",";
        outfile4 << greek.ic_size() << ",";
        outfile4 << greek_anti.mean() << ",";
        outfile4 << greek_anti.ic_size() << ",";
        outfile4 << greek_var_con.mean() << ",";
        outfile4 << greek_var_con.ic_size() << ",";
        outfile4 << greek_var_con_anti.mean() << ",";
        outfile4 << greek_var_con_anti.ic_size() << "\n";

    }
    outfile4.close();
    std::cout << "Results written to Malliavin_BS.csv" << std::endl;
    /*cout<<"greek "<< greek.mean()<<" IC: "<<greek.ic_size()<<endl;
    cout<<"greek anti "<< greek_anti.mean()<<" IC: "<<greek_anti.ic_size()<<endl;
    cout<<"greek con_var "<< greek_var_con.mean()<<" IC: "<<greek_var_con.ic_size()<<endl;
    cout<<"greek con_var_anti "<< greek_var_con_anti.mean()<<" IC: "<<greek_var_con_anti.ic_size()<<endl;*/

//////////////////

    std::cout<<"-------------------------------                   CEV                     -----------------------------------"<<std::endl;

    n_cev=365;
    double h=(static_cast<double>(1)/(n_cev));

    eps=1e-2;
    sample_size_max=1e4;
    std::cout<<"Sample size CEV: "<<sample_size_max<<std::endl;
    pseudo_CEV cev(S,K,r,1, v,T, h);

    diff_finie<pseudo_CEV> df_cev(cev,eps,n_cev,r,T);

    std::cout<<"-------------------------------    Difference finie, pas constant, CEV  -----------------------------------"<<std::endl;
    std::ofstream outfile5;
    outfile4.open("Diff_finie_cst_cev.csv");
    outfile4 << "Sample Size, Greek, IC, Greek Anti, IC Anti, Greek Con Var, IC Con Var\n";
    for (unsigned sample_size = sample_size_init; sample_size <= sample_size_max ; sample_size+=step) {
        auto greek = monte_carlo(&df_cev, &diff_finie<pseudo_CEV>::delta, sample_size, type_payoff, n_cev);
        auto greek_anti = monte_carlo_anti(&df_cev, &diff_finie<pseudo_CEV>::delta, sample_size, type_payoff, n_cev);
        auto greek_var_con = monte_carlo_con_var(&df_cev, &diff_finie<pseudo_CEV>::delta_con_var, sample_size, type_payoff,
                                            n_cev);
        auto greek_var_con_anti = monte_carlo_con_var_anti(&df_cev, &diff_finie<pseudo_CEV>::delta_con_var, sample_size,
                                                      type_payoff, n_cev);

        outfile5 << sample_size << ",";
        outfile5 << greek.mean() << ",";
        outfile5 << greek.ic_size() << ",";
        outfile5 << greek_anti.mean() << ",";
        outfile5 << greek_anti.ic_size() << ",";
        outfile5 << greek_var_con.mean() << ",";
        outfile5 << greek_var_con.ic_size() << ",";
        outfile5 << greek_var_con_anti.mean() << ",";
        outfile5 << greek_var_con_anti.ic_size() << "\n";


        /*cout << "greek " << greek.mean() << " IC: " << greek.ic_size() << endl;
        cout << "greek anti " << greek_anti.mean() << " IC: " << greek_anti.ic_size() << endl;
        cout << "greek con_var " << greek_var_con.mean() << " IC: " << greek_var_con.ic_size() << endl;
        cout << "greek con_var_anti " << greek_var_con_anti.mean() << " IC: " << greek_var_con_anti.ic_size() << endl;*/
    }

    outfile5.close();
    std::cout << "Results written to Diff_finie_cst_cev.csv" << std::endl;



    std::cout<<"-------------------------------    Difference finie, pas decroissant, CEV   -----------------------------------"<<std::endl;


    std::ofstream outfile6;
    outfile6.open("Diff_finie_dec_CEV.csv");
    outfile6 << "Sample Size, Greek, IC, Greek Anti, IC Anti, Greek Con Var, IC Con Var, Greek Anti Con Var, IC Anti Con Var\n";

    for (unsigned sample_size = sample_size_init; sample_size <= sample_size_max; sample_size += step) {
        auto greek_df_d=df_cev.MC_decrease(sample_size,type_payoff,init_eps,n_cev);
        auto greek_df_d_anti=df_cev.MC_decrease_anti(sample_size,type_payoff,init_eps,n_cev);
        auto greek_df_d_con_var=df_cev.MC_decrease_cont_var(r,T,sample_size,type_payoff,init_eps,n_cev);
        auto greek_var_con_anti = monte_carlo_con_var_anti(&df_cev, &diff_finie<pseudo_CEV>::delta_con_var, sample_size,
                                                           type_payoff, n_cev);

        outfile6 << sample_size << ",";
        outfile6 << greek_df_d[2].mean() << ",";
        outfile6 << greek_df_d[0].ic_size() << ",";
        outfile6 << greek_df_d_anti[2].mean() << ",";
        outfile6 << greek_df_d_anti[0].ic_size() << ",";
        outfile6 << greek_df_d_con_var.mean() << ",";
        outfile6 << greek_df_d_con_var.ic_size() << ",";
        outfile6 << greek_var_con_anti.mean() << ",";
        outfile6<< greek_var_con_anti.ic_size() << std::endl;

        /*cout << "greek " << greek_df_d[2].mean() << " IC: " << greek_df_d[0].ic_size() << endl;
        cout << "greek anti " << greek_df_d_anti[2].mean() << " IC: " << greek_df_d_anti[0].ic_size() << endl;
        cout << "greek con_var " << greek_df_d_con_var.mean() << " IC: " << greek_df_d_con_var.ic_size() << endl;
        cout << "greek con_var_anti " << greek_var_con_anti.mean() << " IC: " << greek_var_con_anti.ic_size() << endl;*/
    }

    outfile6.close();
    // Output a message to the console indicating that the CSV file was written
    std::cout << "Results written to Diff_finie_dec_CEV.csv" << std::endl;

    init_eps=1e-3;

    std::ofstream outfile7;
    outfile7.open("Diff_finie_dec_CEV.csv");
    outfile7 << "Sample Size, Greek, IC, Greek Anti, IC Anti, Greek Con Var, IC Con Var, Greek Anti Con Var, IC Anti Con Var\n";

    for (unsigned sample_size = sample_size_init; sample_size <= sample_size_max ; sample_size+=step) {
        auto greek_df_d = df_cev.MC_decrease(sample_size, type_payoff, init_eps, n_cev);
        auto greek_df_d_anti = df_cev.MC_decrease_anti(sample_size, type_payoff, init_eps, n_cev);
        auto greek_df_d_con_var = df_cev.MC_decrease_cont_var(r, T, sample_size, type_payoff, init_eps, n_cev);
        auto greek_var_con_anti = monte_carlo_con_var_anti(&df_cev, &diff_finie<pseudo_CEV>::delta_con_var, sample_size,
                                                           type_payoff, n_cev);

        outfile7 << sample_size << ",";
        outfile7 << greek_df_d[0].mean() << ",";
        outfile7 << greek_df_d[0].ic_size() << ",";
        outfile7 << greek_df_d_anti[0].mean() << ",";
        outfile7 << greek_df_d_anti[0].ic_size() << ",";
        outfile7 << greek_df_d_con_var.mean() << ",";
        outfile7 << greek_df_d_con_var.ic_size() << ",";
        outfile7 << greek_var_con_anti.mean() << ",";
        outfile7<< greek_var_con_anti.ic_size() << std::endl;

        /*cout << "greek " << greek_df_d[0].mean() << " IC: " << greek_df_d[0].ic_size() << endl;
        cout << "greek anti " << greek_df_d_anti[0].mean() << " IC: " << greek_df_d_anti[0].ic_size() << endl;
        cout << "greek con_var " << greek_df_d_con_var.mean() << " IC: " << greek_df_d_con_var.ic_size() << endl;
        cout << "greek con_var_anti " << greek_var_con_anti.mean() << " IC: " << greek_var_con_anti.ic_size() << endl;*/
    }

    outfile7.close();
    std::cout << "Results written to Diff_finie_dec_CEV" << std::endl;


    // Flot tangent 
    tangent_process<pseudo_CEV> tan_cev(r,sample_size_max,T,K,cev,n_cev);

    std::ofstream outfile8;
    outfile8.open("Flot_Tangent_CEV.csv");
    outfile8 << "Sample Size, Greek, IC, Greek Anti, IC Anti, Greek Con Var, IC Con Var, Greek Anti Con Var, IC Anti Con Var\n";

    std::cout<<"-------------------------------            Flot tangent, CEV             -----------------------------------"<<std::endl;

    for (unsigned sample_size = sample_size_init; sample_size <= sample_size_max ; sample_size+=step) {
        auto greek = monte_carlo(&tan_cev, &tangent_process<pseudo_CEV>::delta, sample_size, type_payoff, n_cev);
        auto greek_anti = monte_carlo_anti(&tan_cev, &tangent_process<pseudo_CEV>::delta, sample_size, type_payoff, n_cev);
        auto greek_var_con = monte_carlo_con_var(&tan_cev, &tangent_process<pseudo_CEV>::delta_con_var, sample_size,
                                            type_payoff, n_cev);
        auto greek_var_con_anti = monte_carlo_con_var_anti(&tan_cev, &tangent_process<pseudo_CEV>::delta_con_var,
                                                      sample_size, type_payoff, n_cev);


        // Write the calculated values to the CSV file
        outfile8 << sample_size << ",";
        outfile8 << greek.mean() << ",";
        outfile8 << greek.ic_size() << ",";
        outfile8 << greek_anti.mean() << ",";
        outfile8 << greek_anti.ic_size() << ",";
        outfile8 << greek_var_con.mean() << ",";
        outfile8 << greek_var_con.ic_size() << ",";
        outfile8 << greek_var_con_anti.mean() << ",";
        outfile8 << greek_var_con_anti.ic_size() << "\n";
    }

    outfile8.close();
    std::cout << "Results written to Flot_Tangent_CEV.csv" << std::endl;


    malliavin<pseudo_CEV> mal_cev(r,sample_size_max,T,K,cev,n_cev);

    std::ofstream outfile9;
    outfile9.open("Malliavin_CEV.csv");
    outfile9 << "Sample Size, Greek, IC, Greek Anti, IC Anti, Greek Con Var, IC Con Var\n";

    std::cout<<"-------------------------------             Malliavin, CEV              -----------------------------------"<<std::endl;

    for (unsigned sample_size = sample_size_init; sample_size <= sample_size_max ; sample_size+=step) {
        auto greek = monte_carlo(&mal_cev, &malliavin<pseudo_CEV>::delta, sample_size, type_payoff, n_cev);
        auto greek_anti = monte_carlo_anti(&mal_cev, &malliavin<pseudo_CEV>::delta, sample_size, type_payoff, n_cev);
        auto greek_var_con = monte_carlo_con_var(&mal_cev, &malliavin<pseudo_CEV>::delta_con_var, sample_size, type_payoff,
                                            n_cev);
        auto greek_var_con_anti = monte_carlo_con_var_anti(&mal_cev, &malliavin<pseudo_CEV>::delta_con_var, sample_size,
                                                      type_payoff, n_cev);

        outfile9 << sample_size << ",";
        outfile9 << greek.mean() << ",";
        outfile9 << greek.ic_size() << ",";
        outfile9 << greek_anti.mean() << ",";
        outfile9 << greek_anti.ic_size() << ",";
        outfile9 << greek_var_con.mean() << ",";
        outfile9 << greek_var_con.ic_size() << ",";
        outfile9 << greek_var_con_anti.mean() << ",";
        outfile9 << greek_var_con_anti.ic_size() << "\n";

        /*cout << "greek " << greek.mean() << " IC: " << greek.ic_size() << endl;
        cout << "greek anti " << greek_anti.mean() << " IC: " << greek_anti.ic_size() << endl;
        cout << "greek con_var " << greek_var_con.mean() << " IC: " << greek_var_con.ic_size() << endl;
        cout << "greek con_var_anti " << greek_var_con_anti.mean() << " IC: " << greek_var_con_anti.ic_size() << endl;*/
    }

    outfile9.close();
    std::cout << "Results written to Malliavin_CEV.csv" << std::endl;

}