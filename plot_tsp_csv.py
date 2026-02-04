import csv
from RoutingPractical import plot_tour

def load_tour_from_csv(filename="tsp.csv"):
    """Load tour and coordinates from CSV file exported by C++ code."""
    tour = []
    coords_dict = {}

    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            order = int(row['order'])
            city = int(row['city'])
            x = float(row['x'])
            y = float(row['y'])
            tour.append(city)
            coords_dict[city] = (x, y)

    # Build coords list indexed by city number
    max_city = max(coords_dict.keys())
    coords = [coords_dict[i] for i in range(max_city + 1)]

    return coords, tour

if __name__ == "__main__":
    coords, tour = load_tour_from_csv("tsp.csv")
    plot_tour(coords, tour, title="TSP Tour from C++")
