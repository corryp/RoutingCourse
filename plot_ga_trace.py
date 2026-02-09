import csv
import matplotlib.pyplot as plt

def plot_ga_trace(filename="ga_log.csv"):
    iter = []
    zbest = []
    zmin = []
    zavg = []
    zmax = []

    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            iter.append(int(row['iter']))
            zbest.append(float(row['zbest']))
            zmin.append(float(row['zmin']))
            zavg.append(float(row['zavg']))
            zmax.append(float(row['zmax']))

    fig, ax = plt.subplots(figsize=(10, 5))

    ax.plot(iter, zmax, label='Pop max', alpha=0.4, linewidth=0.5)
    ax.plot(iter, zavg, label='Pop avg', alpha=0.6, linewidth=0.8)
    ax.plot(iter, zmin, label='Pop min', alpha=0.6, linewidth=0.8)
    ax.plot(iter, zbest, label='Best z', linewidth=1.5)
    ax.set_xlabel('Generation')
    ax.set_ylabel('Objective (z)')
    ax.legend()

    plt.title('GA Trace')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    plot_ga_trace()
