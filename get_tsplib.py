import math

def get_tsplib(name):
    """
    Load a TSPLIB instance and return coords and distance matrix.

    Args:
        name: Instance name (e.g., 'berlin52')

    Returns:
        coords: List of (x, y) tuples
        dist: 2D distance matrix
    """
    filename = f"tsplib/{name}.tsp"

    n = 0
    edge_type = "EUC_2D"
    coords = []

    with open(filename, 'r') as f:
        lines = f.readlines()

    i = 0
    # Parse header
    while i < len(lines):
        line = lines[i].strip()

        if line.startswith("DIMENSION"):
            parts = line.replace(":", " ").split()
            n = int(parts[1])
        elif line.startswith("EDGE_WEIGHT_TYPE"):
            parts = line.replace(":", " ").split()
            edge_type = parts[1]
            if edge_type == "EXPLICIT":
                raise ValueError(f"EXPLICIT edge weight type not yet supported: {name}")
        elif line.startswith("NODE_COORD_SECTION"):
            i += 1
            break
        i += 1

    if n == 0:
        raise ValueError(f"Could not find DIMENSION in file: {filename}")

    # Parse coordinates
    coords = [None] * n
    nodes_read = 0
    while i < len(lines) and nodes_read < n:
        line = lines[i].strip()
        if line == "EOF" or line == "":
            i += 1
            continue

        parts = line.split()
        if len(parts) >= 3:
            node_id = int(parts[0])
            x = float(parts[1])
            y = float(parts[2])
            coords[node_id - 1] = (x, y)  # TSPLIB uses 1-based indexing
            nodes_read += 1
        i += 1

    if nodes_read != n:
        raise ValueError(f"Expected {n} nodes but read {nodes_read}")

    # Calculate distance matrix
    dist = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                dist[i][j] = _calc_distance(coords[i], coords[j], edge_type)

    return coords, dist


def _calc_distance(c1, c2, edge_type):
    """Calculate distance between two coordinates based on edge weight type."""
    x1, y1 = c1
    x2, y2 = c2
    dx = x1 - x2
    dy = y1 - y2

    if edge_type == "EUC_2D":
        return math.sqrt(dx * dx + dy * dy)

    elif edge_type == "ATT":
        rij = math.sqrt((dx * dx + dy * dy) / 10.0)
        tij = round(rij)
        return tij + 1 if tij < rij else tij

    elif edge_type == "CEIL_2D":
        return math.ceil(math.sqrt(dx * dx + dy * dy))

    elif edge_type == "GEO":
        PI = 3.141592
        def to_geo(x):
            deg = int(x)
            minute = x - deg
            return PI * (deg + 5.0 * minute / 3.0) / 180.0

        lat1, lon1 = to_geo(x1), to_geo(y1)
        lat2, lon2 = to_geo(x2), to_geo(y2)
        RRR = 6378.388
        q1 = math.cos(lon1 - lon2)
        q2 = math.cos(lat1 - lat2)
        q3 = math.cos(lat1 + lat2)
        return int(RRR * math.acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0)

    else:
        return math.sqrt(dx * dx + dy * dy)


if __name__ == "__main__":
    # Quick test
    coords, dist = get_tsplib("berlin52")
    print(f"Loaded {len(coords)} cities")
    print(f"First 3 coords: {coords[:3]}")
    print(f"Distance 0->1: {dist[0][1]:.2f}")
