import csv
import matplotlib.pyplot as plt

def plot_sa_trace(filename="sa_log.csv"):
    totit = []
    zbest = []
    z = []
    temp = []

    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            totit.append(int(row['totit']))
            zbest.append(float(row['zbest']))
            z.append(float(row['z']))
            temp.append(float(row['temp']))

    fig, ax1 = plt.subplots(figsize=(10, 5))

    ax1.plot(totit, z, label='Incumbent z', alpha=0.5, linewidth=0.5)
    ax1.plot(totit, zbest, label='Best z', linewidth=1.5)
    ax1.set_xlabel('Iteration')
    ax1.set_ylabel('Objective (z)')
    ax1.legend(loc='upper left')

    ax2 = ax1.twinx()
    ax2.plot(totit, temp, label='Temperature', color='red', alpha=0.4, linewidth=0.5)
    ax2.set_ylabel('Temperature')
    ax2.legend(loc='upper right')

    plt.title('SA Trace')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    plot_sa_trace()
