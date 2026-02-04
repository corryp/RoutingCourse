from RoutingPractical import *

#Extract tours from solution
def extract_subtours(x_vars):
    n = len(x_vars)
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
                    if pulp.value(x_vars[current][j]) == 1:
                        next_city = j
                        break
                if next_city is None or visited[next_city]:
                    break
                tour.append(next_city)
                visited[next_city] = True
                current = next_city
            subtours.append(tour)
    return subtours

# Lazy constraints TSP
def tsp_lazy(coords, dist):
    n = len(coords)
    cities = range(n)

    # Initial model with degree constraints only
    prob = pulp.LpProblem("TSP_Lazy", pulp.LpMinimize)
    x = pulp.LpVariable.dicts('x', (cities, cities), cat='Binary')

    # Objective
    prob += pulp.lpSum(dist[i][j] * x[i][j] for i in cities for j in cities if i != j)

    # Degree constraints
    for i in cities:
        prob += pulp.lpSum(x[i][j] for j in cities if j != i) == 1
        prob += pulp.lpSum(x[j][i] for j in cities if j != i) == 1

    # Iterative solve
    start_time = time.time()
    iteration = 0
    while True:
        iteration += 1
        prob.solve(pulp.PULP_CBC_CMD(msg=0))
        subtours = extract_subtours(x)
        print(f"Iteration {iteration}: Found {len(subtours)} subtours")
        if len(subtours) == 1:
            break
        # Add constraints to eliminate subtours
        for tour in subtours:
            if len(tour) < n:
                prob += pulp.lpSum(x[i][j] for i in tour for j in tour if i != j) <= len(tour) - 1

    cpu_time = time.time() - start_time
    report_results(prob, cpu_time, 'lazy')
    return subtours[0], pulp.value(prob.objective), cpu_time


coords, dist = generate_coordinates(15, seed=1)

# Solve MTZ
tour_mtz, obj_mtz, time_mtz = tsp_mtz(coords, dist)

#Solve Lazy
print('\n')
tour_lazy, obj_lazy, time_lazy = tsp_lazy(coords, dist)
plot_tour(coords, tour_lazy, title="Subtour Tour")

