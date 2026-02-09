import csv
import matplotlib.pyplot as plt

def plot_ts_trace(filename="ts_log.csv"):
    iter = []
    zbest = []
    z = []
    tabu = []

    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            iter.append(int(row['iter']))
            zbest.append(float(row['zbest']))
            z.append(float(row['z']))
            tabu.append(int(row['tabu count']))

    fig, ax1 = plt.subplots(figsize=(10, 5))

    ax1.plot(iter, z, label='Current z', alpha=0.5, linewidth=0.5)
    ax1.plot(iter, zbest, label='Best z', linewidth=1.5)
    ax1.set_xlabel('Iteration')
    ax1.set_ylabel('Objective (z)')
    ax1.legend(loc='upper left')

    ax2 = ax1.twinx()
    ax2.plot(iter, tabu, label='Tabu count', color='red', alpha=0.4, linewidth=0.5)
    ax2.set_ylabel('Tabu count')
    ax2.legend(loc='upper right')

    plt.title('TS Trace')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    plot_ts_trace()
