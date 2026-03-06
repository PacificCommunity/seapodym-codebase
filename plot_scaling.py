import numpy as np
import matplotlib.pyplot as plt

num_ranks = np.array([2, 3, 4, 5, 7, 10, 14])
num_workers = num_ranks - 1
time_manager_s = np.array([88.10, 33.19, 23.05, 18.42, 10.27, 9.39, 8.39])

speedup = time_manager_s[0]/time_manager_s
plt.plot(num_workers, speedup)
plt.xlabel('num workers')
plt.ylabel('speedup')
plt.plot(num_workers, num_workers, 'b--')
plt.title('CPU Genoa 1979-01-15 to 2014-12-15')
plt.show()
