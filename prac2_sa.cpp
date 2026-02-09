#include "tsp_mh_config.hpp"
#include "get_tsplib.hpp"

int main() {
    vector<CityPair> v_dat;
    vector<vector<double>> d_dist, d_xy;

    // Load TSP instance (reads name from tsplib_option.csv)
    int n = get_tsplib(v_dat, d_dist, d_xy, "tsplib_option.csv");
    cout << "Loaded TSP instance with " << n << " cities" << endl;

    // Initialise SA components from config file
    TSPnbrOp* nhd = nullptr;
    TSPsoln* x0 = nullptr;
    SActrl hyper_params = initialise_sa(nhd, x0, n, d_dist, "sa_config.csv");

    cout << "Initial solution z = " << x0->md_z << endl;

    // Create solution object for best found
    TSPsoln x_best(n);

    // Set up log file
    ofstream f_log("sa_log.csv");
    Logger log(f_log, false, hyper_params.mi_log_gap);

    // Run simulated annealing
    clock_t t_start = clock();
    SAsummary stats = gx_simulated_annealing(hyper_params, *nhd, *x0, x_best, &log);
    double d_cpu = (clock() - t_start) / (double)CLOCKS_PER_SEC;

    // Output results
    cout << "\n--- SA Results ---" << endl;
    cout << "Initial z = " << stats.md_z0 << endl;
    cout << "Best z = " << stats.md_zbest << endl;
    cout << "T0 = " << stats.md_t0 << endl;
    cout << "Tfinal = " << stats.md_tfinal << endl;
    cout << "Temp steps = " << stats.mi_tstep_ctr << endl;
    cout << "Neighbours evaluated = " << stats.mi_nbr_ctr << endl;
    cout << "Moves made = " << stats.mi_mv_ctr << endl;
    cout << "CPU time = " << d_cpu << "s" << endl;

    // Output best tour to CSV for plotting
    TSPsoln& x_best_tsp = dynamic_cast<TSPsoln&>(x_best);
    g_output_tour_csv(x_best_tsp.mi_x, d_xy);
    cout << "\nBest tour written to tsp.csv" << endl;

    // Cleanup
    delete nhd;
    delete x0;

    return 0;
}
