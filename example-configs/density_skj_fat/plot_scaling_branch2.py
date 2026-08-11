import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import defopt


def main(*, input_file: str='results/timings.csv', output_file: str='results/scaling.png', show: bool=False):
        """
        input_file: CSV input file (e.g. results/timings.csv)
        output_file: PNG output file
        show: set to True to display
        """

        # Load CSV
        df = pd.read_csv(input_file)


        # Ensure manager_time column exists
        if "manager_time" not in df.columns:
            print(f"ERROR: 'manager_time' column not found in CSV. Cannot plot speedup.")

        # Drop rows with missing manager_time
        df = df.dropna(subset=["manager_time"])
        if df.empty:
            print(f"ERROR: No valid 'manager_time' values found. Cannot plot speedup.")
        
        # Plot curves for main and best branches
        cols = {
                 'default': 'r',
                 '-no-aplus-spawn': 'b',
                 '-no-aplus': 'g',
 		}
        plt.figure(figsize=(8,6))
        for host in df["host"].unique():
        
            for branch in df["branch"].unique():
            
                df2 = df[(df.host == host) & (df.branch == branch)]


                df_min = df2.groupby(["nprocs"])["manager_time"].min().reset_index()
                df_avg = df2.groupby(["nprocs"])["manager_time"].mean().reset_index()
                df_max = df2.groupby(["nprocs"])["manager_time"].max().reset_index()
                
                print(f'host = {host} branch = {branch} df_avg = {df_avg}')

                # Compute speedup relative to nprocs=2 for each host, branch 
                df_min['speedup'] = df_min[ df_min["nprocs"] == 2 ]["manager_time"].to_numpy()[0] / df_max["manager_time"].to_numpy()   # note: divide by max
                df_avg['speedup'] = df_avg[ df_min["nprocs"] == 2 ]["manager_time"].to_numpy()[0] / df_avg["manager_time"].to_numpy()
                df_max['speedup'] = df_max[ df_min["nprocs"] == 2 ]["manager_time"].to_numpy()[0] / df_min["manager_time"].to_numpy()   # note: divide by min
            
                col = cols[branch]
                plt.plot(df_min["nprocs"], df_min["speedup"], col + '-.', label=f'{host}-{branch} min')
                plt.plot(df_avg["nprocs"], df_avg["speedup"], col + '-', label=f'{host}-{branch} avg')
                plt.plot(df_max["nprocs"], df_max["speedup"], col + '--', label=f'{host}-{branch} max')

            
        plt.plot([2,20], np.array([2,20]) - 1, 'k-', label='ideal')

        plt.xlabel("Number of processes (nprocs)")
        plt.ylabel("Speedup (relative to nprocs=2)")
        plt.title("Seapodym cohort speedup 1983-2022 skj_fat.xml ")
        plt.xticks(sorted(df["nprocs"].unique()))
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig(output_file)
        if show:
            plt.show()
        print(f"Speedup plot saved to {output_file}")
        
if __name__ == '__main__':
	defopt.run(main)
