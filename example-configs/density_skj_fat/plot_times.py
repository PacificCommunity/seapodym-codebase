import numpy as np
import pandas as pd
import re
import glob
import matplotlib.pyplot as plt
import sys

"""
This script harvests the timing data from the output files results/n*.txt and plots the breakdown of time spent in different components 
(calc, overhead, worker init, cohort init, comm) as a function of the number of ranks. It also computes the average and standard deviation 
of the timings across multiple runs for each component.
"""

pat = re.compile(r'Timings calc/overhead/worker init/cohort init/comm: (\d+\.\d+)\/\s*(\d+\.\d+)\/\s*(\d+\.\d+)\/\s*(\d+\.\d+)\/\s*(\d+\.\d+)')
pat2 = re.compile(r'n(\d+)\.txt')

def get_times(filename: str) -> float:

    times = {
        'calc': [],
        'overhead': [],
        'worker init': [],
        'cohort init': [],
        'comm': []
    }
    f = open (filename, 'r')
    for line in f:
        if 'Timings calc/overhead/worker init/cohort init/comm' in line:
            time_calc, time_overheaad, time_worker_init, time_cohort_init, time_comm = [float(x) for x in re.findall(r'(\d+\.\d+)', line)]
            times['calc'].append(float(time_calc))
            times['overhead'].append(float(time_overheaad))
            times['worker init'].append(float(time_worker_init))
            times['cohort init'].append(float(time_cohort_init))
            times['comm'].append(float(time_comm))
    f.close()
    print(times)
    return times

def main():
    files = glob.glob('results/n*.txt')
    data = {    
        'num_ranks': [],
        'calc': [], 'calc std': [],
        'overhead': [], 'overhead std': [],
        'worker init': [], 'worker init std': [],
        'cohort init': [], 'cohort init std': [],
        'comm': [], 'comm std': [],
    }
    for f in files:
        m = re.search(pat2, f)
        if m:
            n = int(m.group(1))
            data['num_ranks'].append(n)
        else:
            print(f"ERROR Could not extract n from filename {f}")
            sys.exit(1)
        times = get_times(f)
        for name, values in times.items():
            # average the values
            data[name].append(np.mean(values))
            data[name+' std'].append(np.std(values))
    df = pd.DataFrame(data)
    df.sort_values('num_ranks', inplace=True)
    print(df)
    
    df['num workers'] = df['num_ranks'] - 1
    df['ideal'] = 300/df['num workers'] #df[ df['num workers'] == 1 ]['calc'] /df['num workers']
    
    plt.figure()
    for component in ['ideal', 'calc', 'overhead', 'worker init', 'cohort init', 'comm']:
        plt.loglog(df['num workers'], df[component], label=component)
        if hasattr(df, component + ' std'):
            plt.fill_between(df['num workers'], \
                        df[component] - df[component + ' std'], \
                        df[component] + df[component + ' std'], \
                        color='blue', alpha=0.2)

    plt.ylabel('Time (s)')
    plt.title('SEAPODYM timing breakdown (skl_fat.xml na=50 1983-2022)')
    plt.xticks([1, df['num_ranks'].max()], [1, df['num_ranks'].max()])
    plt.legend(title='Component')
    plt.tight_layout()
    plt.xlabel('num workers')
    plt.ylim(bottom=0)
    plt.grid()
    plt.show()

    df['total'] = df['calc'] + df['overhead'] + df['worker init'] + df['cohort init'] + df['comm']

    # speedup
    df['speedup'] = df[ df['num workers'] == 1 ]['total'].values[0] / df['total']
    plt.figure()
    plt.plot(df['num workers'], df['speedup'], marker='o')
    plt.xlabel('num workers')
    plt.ylabel('Speedup')
    plt.title('SEAPODYM speedup (skl_fat.xml na=50 1983-2022)')
    plt.show()

if __name__ == '__main__':
    main()
