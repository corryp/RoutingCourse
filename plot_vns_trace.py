import csv
import matplotlib.pyplot as plt

def plot_vns_trace(filename="vns_log.csv"):
    iter = []
    k = []
    zbest = []
    z = []
    zshake = []
    zls = []

    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            iter.append(int(row['iter']))
            k.append(int(row['k']))
            zbest.append(float(row['zbest']))
            z.append(float(row['z']))
            zshake.append(float(row['zshake']))
            zls.append(float(row['zls']))

    fig, ax1 = plt.subplots(figsize=(10, 5))

    ax1.plot(iter, zshake, label='Shake z', alpha=0.3, linewidth=0.5)
    ax1.plot(iter, zls, label='Local search z', alpha=0.5, linewidth=0.5)
    ax1.plot(iter, z, label='Current z', alpha=0.6, linewidth=0.8)
    ax1.plot(iter, zbest, label='Best z', linewidth=1.5)
    ax1.set_xlabel('Iteration')
    ax1.set_ylabel('Objective (z)')
    ax1.legend(loc='upper left')

    ax2 = ax1.twinx()
    ax2.plot(iter, k, label='Neighbourhood k', color='red', alpha=0.4, linewidth=0.5)
    ax2.set_ylabel('Neighbourhood (k)')
    ax2.legend(loc='upper right')

    plt.title('VNS Trace')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    plot_vns_trace()
