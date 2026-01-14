
"""
TSP with Lazy Constraints using CPLEX Callbacks
------------------------------------------------
This script:
- Generates random coordinates for cities
- Builds TSP model with degree constraints only
- Uses a lazy constraint callback to eliminate subtours dynamically
"""

from docplex.mp.model import Model
import cplex
from cplex.callbacks import LazyConstraintCallback
import time
from RoutingPractical import generate_coordinates, plot_tour

# -----------------------------
# Subtour detection
# -----------------------------
def find_subtours(solution, n):
    visited = [False] * n
    subtours = []
    for start in range(n):
        if not visited[start]:
            tour = [start]
            visited[start] = True
            current = start
            while True:
                next_city = None
                for j in range(n):
                    if solution.get((current, j), 0) > 0.5:
                        next_city = j
                        break
                if next_city is None or visited[next_city]:
                    break
                tour.append(next_city)
                visited[next_city] = True
                current = next_city
            subtours.append(tour)
    return subtours

# -----------------------------
# Extract tour from DOcplex solution
# -----------------------------
def extract_tour_cplex(x, n, mdl):
    tour = [0]
    current = 0
    while True:
        for j in range(n):
            if j != current and mdl.solution.get_value(x[current][j]) > 0.5:
                if j == 0:
                    return tour
                tour.append(j)
                current = j
                break
    return tour

# -----------------------------
# MTZ formulation (CPLEX)
# -----------------------------
def tsp_mtz_cplex(coords, dist):
    n = len(coords)
    mdl = Model("TSP_MTZ")

    # Decision variables
    x = [[mdl.binary_var(name=f"x_{i}_{j}") for j in range(n)] for i in range(n)]
    u = [mdl.continuous_var(lb=0, name=f"u_{i}") for i in range(n)]

    # Objective
    mdl.minimize(mdl.sum(dist[i][j] * x[i][j] for i in range(n) for j in range(n) if i != j))

    # Degree constraints
    for i in range(n):
        mdl.add_constraint(mdl.sum(x[i][j] for j in range(n) if j != i) == 1)
        mdl.add_constraint(mdl.sum(x[j][i] for j in range(n) if j != i) == 1)

    # MTZ subtour elimination constraints
    for i in range(n):
        for j in range(n):
            if i != j and i != 0 and j != 0:
                mdl.add_constraint(u[j] >= u[i] + 1 - (n - 1) * (1 - x[i][j]))
    mdl.add_constraint(u[0] == 0)

    # Solve
    mdl.parameters.timelimit = 60
    start_time = time.time()
    mdl.solve()
    cpu_time = time.time() - start_time

    # Extract tour
    tour = extract_tour_cplex(x, n, mdl)
    obj = mdl.objective_value

    print(f"\n--- MTZ (CPLEX) Results ---")
    print(f"Objective: {obj:.2f}")
    print(f"Variables: {mdl.number_of_variables}")
    print(f"Constraints: {mdl.number_of_constraints}")
    print(f"CPU Time: {cpu_time:.4f} seconds")

    return tour, obj, cpu_time

# -----------------------------
# Lazy constraint callback
# -----------------------------
class SubtourLazyCallback(LazyConstraintCallback):
    def __call__(self):
        n = self.n
        x_indices = self.x_indices

        # Get current solution values
        values = self.get_values(list(x_indices.values()))
        sol = {key: values[idx] for idx, key in enumerate(x_indices.keys())}

        subtours = find_subtours(sol, n)
        if len(subtours) > 1:
            for tour in subtours:
                if len(tour) < n:
                    # Build constraint: sum of x[i][j] for i,j in tour <= |tour| - 1
                    indices = [x_indices[(i, j)] for i in tour for j in tour if i != j]
                    coeffs = [1.0] * len(indices)
                    self.add(constraint=cplex.SparsePair(indices, coeffs),
                            sense="L",
                            rhs=len(tour) - 1)

# -----------------------------
# Build and solve TSP with lazy constraints
# -----------------------------
def tsp_lazy_cplex(coords, dist):
    n = len(coords)
    mdl = Model("TSP_Lazy")

    # Decision variables
    x = [[mdl.binary_var(name=f"x_{i}_{j}") for j in range(n)] for i in range(n)]

    # Objective
    mdl.minimize(mdl.sum(dist[i][j] * x[i][j] for i in range(n) for j in range(n) if i != j))

    # Degree constraints
    for i in range(n):
        mdl.add_constraint(mdl.sum(x[i][j] for j in range(n) if j != i) == 1)
        mdl.add_constraint(mdl.sum(x[j][i] for j in range(n) if j != i) == 1)

    # Export to CPLEX and set up callback
    cpx = mdl.get_cplex()

    # Build mapping from (i,j) to CPLEX variable indices
    x_indices = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                var_name = f"x_{i}_{j}"
                x_indices[(i, j)] = cpx.variables.get_indices(var_name)

    # Register callback with required data
    cb = cpx.register_callback(SubtourLazyCallback)
    cb.n = n
    cb.x_indices = x_indices

    # Solve
    cpx.parameters.timelimit.set(60)
    start_time = time.time()
    cpx.solve()
    cpu_time = time.time() - start_time

    # Extract solution
    solution_values = {}
    for (i, j), idx in x_indices.items():
        solution_values[(i, j)] = cpx.solution.get_values(idx)

    subtours = find_subtours(solution_values, n)
    tour = subtours[0]

    obj = cpx.solution.get_objective_value()

    print(f"\n--- Lazy Constraints (CPLEX) Results ---")
    print(f"Objective: {obj:.2f}")
    print(f"Variables: {mdl.number_of_variables}")
    print(f"Constraints (initial): {mdl.number_of_constraints}")
    print(f"CPU Time: {cpu_time:.4f} seconds")

    return tour, obj, cpu_time

# -----------------------------
# Example usage
# -----------------------------
if __name__ == "__main__":
    coords, dist = generate_coordinates(25, seed=1)

    # Solve with MTZ formulation
    tour_mtz, obj_mtz, time_mtz = tsp_mtz_cplex(coords, dist)

    # Solve with lazy constraints
    tour_lazy, obj_lazy, time_lazy = tsp_lazy_cplex(coords, dist)

    # Summary comparison
    print(f"\n{'='*40}")
    print("COMPARISON SUMMARY")
    print(f"{'='*40}")
    print(f"{'Method':<25} {'Objective':>10} {'Time (s)':>10}")
    print(f"{'-'*40}")
    print(f"{'MTZ':<25} {obj_mtz:>10.2f} {time_mtz:>10.4f}")
    print(f"{'Lazy Constraints':<25} {obj_lazy:>10.2f} {time_lazy:>10.4f}")
    print(f"{'='*40}")

    # Plot the lazy constraints tour
    plot_tour(coords, tour_lazy, title="Lazy Constraints Tour (CPLEX)")
