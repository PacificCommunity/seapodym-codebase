import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import defopt


def main(*, input_file: str='results/timings.csv', output_file: str='results/timings.png', show: bool=False):
        """
        input_file: CSV input file (e.g. results/timings.csv)
        output_file: PNG output file
        show: set to True to display
        """

        # Load CSV
        df = pd.read_csv(input_file)

        # Ensure manager_time column exists
        if "manager_time" not in df.columns:
            print(f"ERROR: 'manager_time' column not found in CSV. Cannot plot.")

        # Drop rows with missing manager_time
        df = df.dropna(subset=["manager_time"])
        if df.empty:
            print(f"ERROR: No valid 'manager_time' values found. Cannot plot.")
        
        # Plot curves for Milan & Genoa
        cols = {
 		'aplus': 'b', 
 		'aplus3': 'm',
 		}
        plt.figure(figsize=(8,6))
        for case in df["case"].unique():
            df2 = df[df.case == case]
            df_min = df2.groupby(["nprocs"])["manager_time"].min().reset_index()
            df_avg = df2.groupby(["nprocs"])["manager_time"].mean().reset_index()
            df_max = df2.groupby(["nprocs"])["manager_time"].max().reset_index()
            
            col = cols[case]
            plt.loglog(df_min["nprocs"], df_min["manager_time"], col + '-.', label=case + ' min')
            plt.loglog(df_avg["nprocs"], df_avg["manager_time"], col + '-', label=case + ' avg')
            plt.loglog(df_max["nprocs"], df_max["manager_time"], col + '--', label=case + ' max')
                        
        plt.loglog([1,100], 1000./np.array([1,100]), 'k-', label='ideal')

        plt.xlabel("Number of processes (nprocs)")
        plt.ylabel("Manager time s")
        plt.title("Seapodym cohort timings 1983-2022 skj_Fat.xml A+")
        plt.xticks(sorted(df["nprocs"].unique()))
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig(output_file)
        if show:
            plt.show()
        print(f"plot saved to {output_file}")
        
if __name__ == '__main__':
	defopt.run(main)
