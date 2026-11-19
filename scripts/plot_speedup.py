
import numpy as np
import matplotlib.pyplot as plt

# MPI process counts
processes = np.array([1, 2, 4, 16])

# Total runtime in seconds, 1000 iterations per trial
jacobi_trials = [
    [0.00725075, 0.00721433, 0.00740147, 0.00723762, 0.00710558],
    [0.00503716, 0.00489926, 0.00489285, 0.00494247, 0.00498266],
    [0.00456627, 0.00451601, 0.00450630, 0.00451970, 0.00455456],
    [0.01050250, 0.00761411, 0.00752560, 0.02894430, 0.00726936],
]

rbgs_trials = [
    [0.0318378, 0.0318939, 0.0317462, 0.0320292, 0.0319254],
    [0.0187463, 0.0184510, 0.0330719, 0.0183328, 0.0186858],
    [0.0139917, 0.0132336, 0.0132705, 0.0135737, 0.0133411],
    [0.0212607, 0.0175896, 0.0178061, 0.0479263, 0.0342006],
]

# Median runtime for each process count
jacobi_time = np.median(jacobi_trials, axis=1)
rbgs_time = np.median(rbgs_trials, axis=1)

# Strong scaling speedup
jacobi_speedup = jacobi_time[0] / jacobi_time
rbgs_speedup = rbgs_time[0] / rbgs_time

plt.figure(figsize=(8, 5))

plt.plot(processes, processes, "k--", label="Ideal linear speedup")
plt.plot(processes, jacobi_speedup, "o-", label="Jacobi")
plt.plot(processes, rbgs_speedup, "s-", label="Red-black Gauss-Seidel")

plt.xlabel("Number of MPI processes")
plt.ylabel("Speedup relative to 1 process")
plt.title("Strong Scaling of MPI Iterative Solvers (128 x 128)")
plt.xticks(processes)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig("speedup.png", dpi=300)
plt.show()

print("MPI ranks | Jacobi speedup | RBGS speedup")
for p, sj, sr in zip(processes, jacobi_speedup, rbgs_speedup):
    print(f"{p:9d} | {sj:14.3f} | {sr:12.3f}")
