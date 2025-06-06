
import matplotlib.pyplot as plt
import csv
import os
import numpy as np

with open("positions_sequential.csv") as f:
    reader = csv.reader(f)

    mass_row = next(reader)
    masses = [float(val) for val in mass_row if val.strip() != '']

    data = []
    for row in reader:
        cleaned_row = [val for val in row if val.strip() != '']
        if cleaned_row:
            try:
                data.append(list(map(float, cleaned_row)))
            except ValueError:
                continue

data = np.array(data)
num_bodies = len(masses)

mass_array = np.array(masses)
sizes = 10 + 40 * np.log10(mass_array / mass_array.min())

colors = ['yellow'] + [plt.cm.tab10(i) for i in range(num_bodies - 1)] # yellow sun


os.makedirs("frames_mass", exist_ok=True)

for i in range(data.shape[0]):
    fig, ax = plt.subplots(figsize=(6, 6), facecolor='black')
    ax.set_facecolor('black')

    x_vals = [data[i, 2*b] for b in range(num_bodies)]
    y_vals = [data[i, 2*b+1] for b in range(num_bodies)]

    x_min, x_max = min(x_vals), max(x_vals)
    y_min, y_max = min(y_vals), max(y_vals)
    margin = 0.05 * max(x_max - x_min, y_max - y_min, 1e10)

    ax.set_xlim(x_min - margin, x_max + margin)
    ax.set_ylim(y_min - margin, y_max + margin)

    for b in range(num_bodies):
        x = data[i, 2 * b]
        y = data[i, 2 * b + 1]
        ax.scatter(x, y, s=sizes[b], color=colors[b], label=f"Body {b}" if i == 0 else "")

    if i == 0:
        ax.legend(loc="upper right", fontsize=8)

    ax.set_title(f"Step {i}", color='white')
    ax.tick_params(colors='white')  
    for spine in ax.spines.values():
        spine.set_edgecolor('white')

    plt.savefig(f"frames_mass/frame_{i:03d}.png", facecolor=fig.get_facecolor())
    plt.close()

