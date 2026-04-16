import re
import glob
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import defopt


COLORS = {
    "init": "green",
    "step": "lightblue",
    "put": "red",
    "notify": "blue",
}

# ---------------- PARSER ----------------
def parse_logs(pattern):

    files = sorted(glob.glob(pattern))
    print(f"Found {len(files)} files")

    records = []
    active = {}

    ts_re = r"\[(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d+)\]"
    worker_re = r"\[log(\d+)\]"
    task_re = r"task id (\d+)"
    step_re = r"step (\d+)"

    task_ids = []
    worker_ids = []
    phases = []
    t_starts = []
    t_ends = []
    steps = []

    for file in files:
        print(f"Parsing {file}")

        with open(file) as f:

            phase = None
            t_start = None

            for line in f:

                # ---------- timestamp ----------
                ts_m = re.search(ts_re, line)
                if not ts_m:
                    continue

                ts = datetime.strptime(ts_m.group(1), "%Y-%m-%d %H:%M:%S.%f")
                worker_id = int(re.search(worker_re, line).group(1))

                if phase:

                    m = re.search(r'<< initialization of task id (\d+)', line)
                    if m:
                        
                        # end of init
                        task_id = int(m.group(1))
                        task_ids.append(task_id)
                        worker_ids.append(worker_id)
                        phases.append('init')
                        t_starts.append(t_start)
                        t_ends.append(ts)
                        steps.append(-1)
                        phase = None
                        continue

                    m = re.search(r'<<< send data for step (\d+) of task id (\d+)', line)
                    if m:
                        # end of put
                        task_id = int(m.group(2))
                        task_ids.append(task_id)
                        worker_ids.append(worker_id)
                        phases.append('put')
                        t_starts.append(t_start)
                        t_ends.append(ts)
                        steps.append(-2)
                        phase = None
                        continue

                    m = re.search(r'<<< step (\d+) of task id (\d+)', line)
                    if m:
                        # end of step
                        task_id = int(m.group(2))
                        task_ids.append(task_id)
                        worker_ids.append(worker_id)
                        phases.append('step')
                        t_starts.append(t_start)
                        t_ends.append(ts)
                        steps.append(int(m.group(1)))
                        phase = None
                        continue

                    m = re.search(r'<<< notify manager after step (\d+) of task id (\d+)', line)
                    if m:
                        # end of notify
                        task_id = int(m.group(2))
                        task_ids.append(task_id)
                        worker_ids.append(worker_id)
                        phases.append('notify')
                        t_starts.append(t_start)
                        t_ends.append(ts)
                        steps.append(-3)
                        phase = None
                        continue

                else:
                
                    m = re.search(r'>> initialization of task id (\d+)', line)
                    if m:
                        t_start = ts
                        phase = 'init'
                        continue

                    m = re.search(r'>>> send data', line)
                    if m:
                        t_start = ts
                        phase = 'put'
                        continue

                    m = re.search(r'>>> step (\d+) of task id (\d+)', line)
                    if m:
                        t_start = ts
                        phase = 'step'
                        continue

                    m = re.search(r'>>> notify manager after step (\d+) of task id (\d+)', line)
                    if m:
                        t_start = ts
                        phase = 'notify'
                        continue
 
    df = pd.DataFrame({
        'task_id': task_ids,
        'worker_id': worker_ids,
        'phase': phases,
        'step': steps,
        't_start': t_starts,
        't_end': t_ends,
    })

    print("\nPhases found:", df["phase"].unique())
    print("Total rows:", len(df))
    print(df)

    return df


# ---------------- PLOT ----------------
def plot_gantt(df):
    fig, ax = plt.subplots(figsize=(14, 6))

    for _, r in df.iterrows():
        ax.barh(
            r["worker_id"],
            r["end_s"] - r["start_s"],
            left=r["start_s"],
            height=1.0, # fill vertical space, no gaps
            color=COLORS.get(r["phase"], "gray"),
        )

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Worker ID")
    ax.set_title("Worker Timeline")

    ax.set_yticks(sorted(df["worker_id"].unique()))

    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in COLORS.values()]
    ax.legend(handles, list(COLORS.keys()))

    plt.tight_layout()
    plt.show()

def print_summary(df):
    total_time = df.end_s.max()
    num_workers = df.worker_id.max()
    for worker_id in range(1, num_workers):
        df2 = df[ df['worker_id'] == worker_id ]
        print(f'\nWorker {worker_id}')
        for phase in df.phase.unique():
            df3 = df2[ df2['phase'] == phase]
            dt = df3['end_s'] - df3['start_s']
            dt_sum = dt.sum()
            total_time_worker = df3.end_s.max() - df3.start_s.min()
            print(f'{phase}\t: {dt.min():.3f} <= {dt.mean():.3f} +/- {dt.std():.3f} <= {dt.max():.3f} total={dt_sum:.3f} {dt_sum*100/(total_time_worker):.1f}%')
    print('\nAll workers')
    for phase in df.phase.unique():
        df3 = df[ df['phase'] == phase]
        dt = df3['end_s'] - df3['start_s']
        dt_sum = dt.sum()
        total_time = df3.end_s.max() - df3.start_s.min()
        print(f'{phase}\t: {dt.min():.3f} <= {dt.mean():.3f} +/- {dt.std():.3f} <= {dt.max():.3f} total={dt_sum:.3f} {dt_sum*100/(num_workers*total_time):.1f}%')
     

# ---------------- MAIN ----------------
def main(*, log_pattern: str="log_taskfunc*.txt", tmin: float=0, tmax: float=-1, show: bool=False):

    df = parse_logs(log_pattern)

    if not df.empty:
        df = df.sort_values(["worker_id", "t_start"])

    t0 = df["t_start"].min()
    df["start_s"] = (df["t_start"] - t0).dt.total_seconds()
    df["end_s"] = (df["t_end"] - t0).dt.total_seconds()

    exec_time = df["end_s"].max()
    print(f'exec time: {exec_time}')
    if tmax < 0:
        tmax = exec_time
    
    # select time interval
    df = df[df.start_s >= tmin]
    df = df[df.end_s < tmax]

    print_summary(df)

    if show:
    	plot_gantt(df)

if __name__ == "__main__":
    defopt.run(main)

