

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


colors = plt.cm.get_cmap('tab10', num_bodies)


os.makedirs("frames_mass", exist_ok=True)

for i in range(data.shape[0]):
    plt.figure(figsize=(6,6))

    x_vals = [data[i, 2*b] for b in range(num_bodies)]
    y_vals = [data[i, 2*b+1] for b in range(num_bodies)]

    x_min, x_max = min(x_vals), max(x_vals)
    y_min, y_max = min(y_vals), max(y_vals)
    margin = 0.1 * max(x_max - x_min, y_max - y_min, 1e10)

    plt.xlim(x_min - margin, x_max + margin)
    plt.ylim(y_min - margin, y_max + margin)

    for b in range(num_bodies):
        x = data[i, 2 * b]
        y = data[i, 2 * b + 1]
        plt.scatter(x, y, s=sizes[b], color=colors(b), label=f"Body {b}" if i == 0 else "")

    plt.title(f"Step {i}")
    if i == 0:
        plt.legend(loc="upper right")
    plt.savefig(f"frames_mass/frame_{i:03d}.png")
    plt.close()


