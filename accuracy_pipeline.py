import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import os


if os.path.exists("accuracy_results.csv"):
    os.remove("accuracy_results.csv")

with open("accuracy_results.csv", "w") as f:
    f.write("num_bodies,parallel,mutex,nomutex,barneshutt\n")

bodiesList = [50, 100, 1000, 2000, 5000, 10000]

for bodies in bodiesList:
    print(f"Running simulation with {bodies} bodies...")
    subprocess.run(["./simple_approach_test", str(bodies)])


df = pd.read_csv("accuracy_results.csv")

plt.figure(figsize=(10, 6))
plt.plot(df['num_bodies'], df['parallel'], 'o-', label='Parallel')
plt.plot(df['num_bodies'], df['mutex'], 'o-', label='Parallel + Mutex')
plt.plot(df['num_bodies'], df['nomutex'], 'o-', label='Parallel No Mutex')
plt.plot(df['num_bodies'], df['barneshutt'], 'o-', label='Barnes-Hutt')

plt.xlabel("Number of Bodies")
plt.ylabel("Accuracy (%) vs. Sequential")
plt.title("Accuracy of Methods vs. Number of Bodies, theta = 0.9")
plt.ylim(99.5, 100.0)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("accuracy_plot.png")
plt.show()
