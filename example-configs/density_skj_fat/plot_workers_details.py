import re
import glob
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

LOG_PATTERN = "log_taskfunc1.txt"

COLORS = {
    "init": "green",
    "step": "lightblue",
    "put": "red",
    "notify": "orange",
}

re_init_start = re.compile(r"")

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

                print(f'**** line = {line}')

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
                        print(f'<< init detected at line: {line} ts = {ts}')
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
                        print(f'<<< put detected at line: {line} ts = {ts}')
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
                        print(f'>> init detected at line: {line} ts = {ts}')
                        continue

                    m = re.search(r'>>> send data', line)
                    if m:
                        t_start = ts
                        phase = 'put'
                        print(f'>>> put detected at line: {line}')
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

    #print("\nPhases found:", df["phase"].unique())
    print("Total rows:", len(df))
    print(df)

    return df


# ---------------- PLOT ----------------
def plot_gantt(df):
    fig, ax = plt.subplots(figsize=(14, 6))

    t0 = df["t_start"].min()
    df["start_s"] = (df["t_start"] - t0).dt.total_seconds()
    df["end_s"] = (df["t_end"] - t0).dt.total_seconds()

    for _, r in df.iterrows():
        ax.barh(
            r["worker_id"],
            r["end_s"] - r["start_s"],
            left=r["start_s"],
            color=COLORS.get(r["phase"], "gray"),
        )

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Worker ID")
    ax.set_title("MPI Worker Timeline (Corrected Parser)")

    ax.set_yticks(sorted(df["worker_id"].unique()))

    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in COLORS.values()]
    ax.legend(handles, list(COLORS.keys()))

    plt.tight_layout()
    plt.show()


# ---------------- MAIN ----------------
if __name__ == "__main__":
    df = parse_logs(LOG_PATTERN)

    print(df.head(20))
    print(df.phase.unique())

    # if not df.empty:
    #     df = df.sort_values(["worker_id", "t_start"])

    # plot_gantt(df)



