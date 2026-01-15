#include "BasicTSPheuristics.hpp"
#include "get_tsplib.hpp"

void solve_tsplib(string tsp_name) {
	vector<CityPair> v_dat;
	vector<vector<double>> d_dist, d_xy;
    
    int n = get_tsplib(tsp_name, v_dat, d_dist, d_xy);
    vector<int> i_rnd_tour(n);
    double d_zrnd = gd_random_tour(n, d_dist, i_rnd_tour);

    vector<int> i_nn_tour(n);
    double d_znn = gd_nearest_neighbour(n, d_dist, 0, i_nn_tour);

    vector<int> i_2o_rnd_tour(i_rnd_tour);
    vector<int> i_2o_nn_tour(i_nn_tour);

    clock_t t_start_rnd = clock();
    double d_z2o_rnd = gd_local_2opt_search(n, d_dist, i_2o_rnd_tour, true);
    double d_cpu_rnd = (clock() - t_start_rnd) / (double)CLOCKS_PER_SEC;

    clock_t t_start_nn = clock();
    double d_z2o_nn = gd_local_2opt_search(n, d_dist, i_2o_nn_tour, true);
    double d_cpu_nn = (clock() - t_start_nn) / (double)CLOCKS_PER_SEC;

    cout << "z_rnd = " << d_zrnd << endl;
    cout << "z_nn = " << d_znn << endl;
    cout << "z_2o_rnd = " << d_z2o_rnd << " (CPU: " << d_cpu_rnd << "s)" << endl;
    cout << "z_2o_nn = " << d_z2o_nn << " (CPU: " << d_cpu_nn << "s)" << endl;

    g_output_tour_csv(i_2o_nn_tour, d_xy);
}

int main() {
    solve_tsplib("ch130");
    return 0;
}