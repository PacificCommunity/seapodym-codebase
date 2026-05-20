import numpy as np
import pandas as pd
import re
import glob
import matplotlib.pyplot as plt
import sys

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
            time_calc, time_overheaad, time_worker_init, time_cohort_init, time_comm = [float(x) for x in re.findall('(\d+\.\d+)', line)]
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
    
    plt.figure()
    for component in ['overhead', 'worker init', 'cohort init', 'comm']:
        plt.plot(df['num_ranks'], df[component], label=component)
    plt.plot(df['num_ranks'], df['calc'], label='calc')
    plt.ylabel('Time (s)')
    plt.title('Timing breakdown by component')
    plt.xticks(rotation=0)
    plt.legend(title='Component')
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()