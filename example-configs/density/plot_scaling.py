import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


if __name__ == '__main__':

        # Load CSV
        df = pd.read_csv('results/timings.csv')


        # Ensure manager_time column exists
        if "manager_time" not in df.columns:
            print(f"ERROR: 'manager_time' column not found in CSV. Cannot plot speedup.")

        # Drop rows with missing manager_time
        df = df.dropna(subset=["manager_time"])
        if df.empty:
            print(f"ERROR: No valid 'manager_time' values found. Cannot plot speedup.")

        # Take minimum manager_time per (nprocs, host)
        df_min = df.groupby(["nprocs","host"])["manager_time"].min().reset_index()

        # Compute speedup relative to nprocs=2 for each host
        df_speedup_list = []
        for host in df_min["host"].unique():
            host_df = df_min[df_min["host"] == host].copy()
            if not 2 in host_df["nprocs"].values:
                print(f"WARNING: nprocs=2 missing for host '{host}', skipping speedup calculation for this host.")
                continue
            t2 = host_df.loc[host_df["nprocs"]==2, "manager_time"].values[0]
            host_df["speedup"] = t2 / host_df["manager_time"]
            df_speedup_list.append(host_df)

        if not df_speedup_list:
            print("ERROR: No valid hosts with nprocs=2. Cannot plot speedup.")

        df_speedup = pd.concat(df_speedup_list)

        # Plot curves for Milan & Genoa
        plt.figure(figsize=(8,6))
        for host in ["milan","genoa"]:
            host_df = df_speedup[df_speedup["host"]==host]
            if not host_df.empty:
                plt.plot(host_df["nprocs"], host_df["speedup"], marker='o', label=host.capitalize())
        plt.plot([2,20], np.array([2,20]) - 1, 'k--', label='ideal')

        plt.xlabel("Number of processes (nprocs)")
        plt.ylabel("Speedup (relative to nprocs=2)")
        plt.title("Seapodym cohort speedup")
        plt.xticks(sorted(df_speedup["nprocs"].unique()))
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig('scaling.png')
        plt.show()
        print(f"Speedup plot saved to scaling.png")
