import pulp
import itertools
import random
import math
import time
import matplotlib.pyplot as plt

# -----------------------------
# Generate coordinates and distance matrix
# -----------------------------
def generate_coordinates(n, seed=None):
    if seed is not None:
        random.seed(seed)
    coords = [(random.uniform(0, 100), random.uniform(0, 100)) for _ in range(n)]
    dist = [[0 if i == j else math.dist(coords[i], coords[j]) for j in range(n)] for i in range(n)]
    return coords, dist

# -----------------------------
# Extract tour from decision variables
# -----------------------------
def extract_tour(x_vars):
    tour = [0]  # start at city 0
    current = 0
    while True:
        next_city = None
        for j in x_vars[current]:
            if pulp.value(x_vars[current][j]) == 1:
                next_city = j
                break
        if next_city == 0:
            break
        tour.append(next_city)
        current = next_city
    return tour

# -----------------------------
# Plot tour
# -----------------------------
def plot_tour(coords, tour, title="TSP Tour"):
    x_vals = [coords[i][0] for i in tour] + [coords[tour[0]][0]]
    y_vals = [coords[i][1] for i in tour] + [coords[tour[0]][1]]
    length = sum(math.dist(coords[tour[i]], coords[tour[(i + 1) % len(tour)]]) for i in range(len(tour)))
    plt.figure(figsize=(6, 6))
    plt.scatter([c[0] for c in coords], [c[1] for c in coords], color='blue')
    plt.plot(x_vals, y_vals, color='red', linewidth=2)
    for idx, (x, y) in enumerate(coords):
        plt.text(x, y, str(idx), fontsize=9)
    plt.title(f"{title} (length: {length:.2f})")
    plt.show()


# -----------------------------
# Report solver results
# -----------------------------
def report_results(prob, cpu_time, formulation_name):
    print(f"\n--- {formulation_name} Results ---")
    print("Status:", pulp.LpStatus[prob.status])
    print("Optimal Objective:", pulp.value(prob.objective))
    print("Variables:", len(prob.variables()))
    print("Constraints:", len(prob.constraints))
    print(f"CPU Time: {cpu_time:.4f} seconds")
    return cpu_time


# -----------------------------
# MTZ formulation
# -----------------------------
def tsp_mtz(coords, dist):
    n = len(coords)
    cities = range(n)
    prob = pulp.LpProblem("TSP_MTZ", pulp.LpMinimize)
    x = pulp.LpVariable.dicts('x', (cities, cities), cat='Binary')
    u = pulp.LpVariable.dicts('u', cities, lowBound=0, cat='Continuous')

    prob += pulp.lpSum(dist[i][j] * x[i][j] for i in cities for j in cities if i != j)

    for i in cities:
        prob += pulp.lpSum(x[i][j] for j in cities if j != i) == 1
        prob += pulp.lpSum(x[j][i] for j in cities if j != i) == 1

    for i in cities:
        for j in cities:
            if i != j and i != 0 and j != 0:
                prob += u[j] >= u[i] + 1 - (n - 1) * (1 - x[i][j])
    prob += u[0] == 0

    start_time = time.time()
    prob.solve(pulp.PULP_CBC_CMD(msg=0))
    cpu_time = time.time() - start_time
    report_results(prob, cpu_time, 'MTZ')

    tour = extract_tour(x)
    return tour, pulp.value(prob.objective), cpu_time

# -----------------------------
# Subtour formulation
# -----------------------------
def tsp_subtour(coords, dist):
    n = len(coords)
    cities = range(n)
    prob = pulp.LpProblem("TSP_Subtour", pulp.LpMinimize)
    x = pulp.LpVariable.dicts('x', (cities, cities), cat='Binary')

    prob += pulp.lpSum(dist[i][j] * x[i][j] for i in cities for j in cities if i != j)

    for i in cities:
        prob += pulp.lpSum(x[i][j] for j in cities if j != i) == 1
        prob += pulp.lpSum(x[j][i] for j in cities if j != i) == 1

    subsets = []
    for r in range(2, n):
        subsets += itertools.combinations(cities, r)
    for S in subsets:
        prob += pulp.lpSum(x[i][j] for i in S for j in S if i != j) <= len(S) - 1

    start_time = time.time()
    prob.solve(pulp.PULP_CBC_CMD(msg=0))
    cpu_time = time.time() - start_time

    report_results(prob, cpu_time, 'subtour constraints')

    tour = extract_tour(x)
    return tour, pulp.value(prob.objective), cpu_time

# -----------------------------
# Example usage
# -----------------------------

if __name__ == "__main__":
    coords, dist = generate_coordinates(10, seed=1)

    # Solve MTZ
    tour_mtz, obj_mtz, time_mtz = tsp_mtz(coords, dist)

    # Solve Subtour
    tour_sub, obj_sub, time_sub = tsp_subtour(coords, dist)
    plot_tour(coords, tour_sub, title="Subtour Tour")



