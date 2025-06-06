
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

planet_names = ["Sun", "Mercury", "Venus", "Earth", "Mars"]

sizes = 10 + 40 * np.log10(mass_array / mass_array[1]) # to normalize sizes
sizes[0] = 30  # smaller Sun

colors = ['yellow'] + [plt.cm.tab10(i) for i in range(num_bodies - 1)]

os.makedirs("frames_mass", exist_ok=True)


focus_bodies = list(range(5)) #zoom on first 4 planets


trails = [[] for _ in focus_bodies]

for i in range(data.shape[0]):
    fig, ax = plt.subplots(figsize=(6, 6), facecolor='black')
    ax.set_facecolor('black')

    sun_x = data[i, 0]
    sun_y = data[i, 1]

    view_radius = 2.5e11 # sun in center
    ax.set_xlim(sun_x - view_radius, sun_x + view_radius)
    ax.set_ylim(sun_y - view_radius, sun_y + view_radius)

    for idx, b in enumerate(focus_bodies):
        x = data[i, 2 * b]
        y = data[i, 2 * b + 1]

        trails[idx].append((x, y))

        trail_array = np.array(trails[idx])
        if len(trail_array) > 1:
            ax.plot(trail_array[:, 0], trail_array[:, 1], color='white', alpha=0.15, linewidth=0.5)


        ax.scatter(x, y, s=sizes[b], color=colors[b])


        ax.text(x + 0.01 * view_radius, y + 0.01 * view_radius, # for the names of planets
                planet_names[idx], fontsize=6, color='white')

    ax.set_title(f"Step {i}", color='white')
    ax.tick_params(colors='white')
    for spine in ax.spines.values():
        spine.set_edgecolor('white')

    plt.savefig(f"frames_mass/frame_{i:03d}.png", facecolor=fig.get_facecolor())
    plt.close()
