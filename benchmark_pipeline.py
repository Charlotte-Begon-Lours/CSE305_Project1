import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import os

if os.path.exists("timing_results.csv"):
    os.remove("timing_results.csv")

with open("timing_results.csv", "w") as f:
    f.write("threads,sequential\n") #,Barnes-Hut,Barnes-Hut partially parallel,Barnes-Hut parallel

for threads in range(1, 17):
    print(f"Running simulation with {threads} thread(s)...")
    subprocess.run(["./simple_approach_test", str(threads)])


df = pd.read_csv("timing_results.csv")

plt.figure(figsize=(10, 6))
plt.plot(df['threads'], df['sequential'], 'o-', label='Sequential')
# plt.plot(df['threads'], df['parallel'], 'o-', label='Parallel')
# plt.plot(df['threads'], df['mutex'], 'o-', label='Parallel + Mutex')
# plt.plot(df['threads'], df['nomutex'], 'o-', label='Parallel No Mutex')
#plt.plot(df['threads'], df['Barnes-Hut'], 'o-', label='Barnes-Hut')
#plt.plot(df['threads'], df['Barnes-Hut partially parallel'], 'o-', label='Barnes-Hut partially parallelized')
#plt.plot(df['threads'], df['Barnes-Hut parallel'], 'o-', label='Parallel Barnes-Hut')

plt.xlabel("Number of Threads")
plt.ylabel("Runtime (ms)")
plt.title("Runtime vs. Thread Count for 50 bodies, theta = 0.5")
plt.legend()
plt.grid(True)
plt.savefig("runtime_plot.png")
plt.show()
