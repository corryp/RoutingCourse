#include "tsp_mh_config.hpp"
#include "get_tsplib.hpp"

int main() {
    vector<CityPair> v_dat;
    vector<vector<double>> d_dist, d_xy;

    // Load TSP instance (reads name from tsplib_option.csv)
    int n = get_tsplib(v_dat, d_dist, d_xy, "tsplib_option.csv");
    cout << "Loaded TSP instance with " << n << " cities" << endl;

    // Initialise GA components from config file
    TSPnbrOp* mutt = nullptr;
    TSPxover* xover = nullptr;
    GActrl hyper_params = initialise_ga(mutt, xover, n, d_dist, "ga_config.csv");

    // Create solution object for best found
    TSPsoln x_best(n);

    // Run genetic algorithm
    clock_t t_start = clock();
    GAtrace stats = gx_genetic_alg(hyper_params, *mutt, *xover, x_best);
    double d_cpu = (clock() - t_start) / (double)CLOCKS_PER_SEC;

    // Output results
    cout << "\n--- GA Results ---" << endl;
    cout << "Best z = " << stats.md_zbest << endl;
    cout << "Final pop min z = " << stats.md_zmin << endl;
    cout << "Final pop avg z = " << stats.md_zavg << endl;
    cout << "Final pop max z = " << stats.md_zmax << endl;
    cout << "Generations = " << stats.mi_gen_ctr << endl;
    cout << "CPU time = " << d_cpu << "s" << endl;

    // Output best tour to CSV for plotting
    TSPsoln& x_best_tsp = dynamic_cast<TSPsoln&>(x_best);
    g_output_tour_csv(x_best_tsp.mi_x, d_xy);
    cout << "\nBest tour written to tsp.csv" << endl;

    // Cleanup
    delete mutt;
    delete xover;

    return 0;
}
