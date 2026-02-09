#include "tsp_mh_config.hpp"
#include "get_tsplib.hpp"

int main() {
    vector<CityPair> v_dat;
    vector<vector<double>> d_dist, d_xy;

    // Load TSP instance (reads name from tsplib_option.csv)
    int n = get_tsplib(v_dat, d_dist, d_xy, "tsplib_option.csv");
    cout << "Loaded TSP instance with " << n << " cities" << endl;

    // Initialise TS components from config file
    TSPnbrOp* nhd = nullptr;
    TSPsoln* x0 = nullptr;
    TabuListZhash* tabu_list = nullptr;
    TSPnhdMem* long_mem = nullptr;
    TSctrl hyper_params = initialise_ts(nhd, x0, tabu_list, long_mem, n, d_dist, "ts_config.csv");

    cout << "Initial solution z = " << x0->md_z << endl;

    // Create solution object for best found
    TSPsoln x_best(n);

    // Set up log file
    ofstream f_log("ts_log.csv");
    Logger log(f_log, false, hyper_params.mi_log_gap);

    // Run tabu search
    clock_t t_start = clock();
    TSsummary stats = gx_tabu_search(hyper_params, *tabu_list, long_mem, *nhd, *x0, x_best, &log);
    double d_cpu = (clock() - t_start) / (double)CLOCKS_PER_SEC;

    // Output results
    cout << "\n--- TS Results ---" << endl;
    cout << "Initial z = " << stats.md_z0 << endl;
    cout << "Best z = " << stats.md_zbest << endl;
    cout << "Iterations = " << stats.mi_nhd_ctr << endl;
    cout << "CPU time = " << d_cpu << "s" << endl;

    // Output best tour to CSV for plotting
    TSPsoln& x_best_tsp = dynamic_cast<TSPsoln&>(x_best);
    g_output_tour_csv(x_best_tsp.mi_x, d_xy);
    cout << "\nBest tour written to tsp.csv" << endl;

    // Cleanup
    delete nhd;
    delete x0;
    delete tabu_list;
    delete long_mem;

    return 0;
}
