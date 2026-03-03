import pandas as pd
import re
from datetime import datetime
from glob import glob
import matplotlib.pyplot as plt


def parse_logs():

    # --- Regular expressions ---
    re_init_start = re.compile(r"\[(.*?)\].*>> initialization of task id (\d+)")
    re_init_end   = re.compile(r"\[(.*?)\].*<< initialization of task id (\d+)")
    re_step_start = re.compile(r"\[(.*?)\].*>>> step (\d+) of task id (\d+)")
    re_step_end   = re.compile(r"\[(.*?)\].*<<< step (\d+) of task id (\d+)")
    re_send_start = re.compile(r"\[(.*?)\].*>>> send data for step (\d+) of task id (\d+)")
    re_send_end   = re.compile(r"\[(.*?)\].*<<< send data for step (\d+) of task id (\d+)")
    re_worker_id  = re.compile(r"log_worker(\d+).txt")

    # --- Helper to convert timestamps ---
    def parse_time(s):
        # e.g. "2025-10-22 03:52:18.536"
        return datetime.strptime(s, "%Y-%m-%d %H:%M:%S.%f")

    # --- Collect raw data ---
    records = []

    for fname in sorted(glob("log_worker*.txt")):

        worker_id = None
        if m := re_worker_id.search(fname):
            worker_id = int(m[1])

        with open(fname) as f:
            lines = f.readlines()

        init_start, init_end = {}, {}
        step_start, step_end = {}, {}
        send_start, send_end = {}, {}

        for line in lines:
            if m := re_init_start.search(line):
                t, task = parse_time(m[1]), int(m[2])
                init_start[task] = t
            elif m := re_init_end.search(line):
                t, task = parse_time(m[1]), int(m[2])
                init_end[task] = t
            elif m := re_step_start.search(line):
                t, step, task = parse_time(m[1]), int(m[2]), int(m[3])
                step_start[(task, step)] = t
            elif m := re_step_end.search(line):
                t, step, task = parse_time(m[1]), int(m[2]), int(m[3])
                step_end[(task, step)] = t
            elif m := re_send_start.search(line):
                t, step, task = parse_time(m[1]), int(m[2]), int(m[3])
                send_start[(task, step)] = t
            elif m := re_send_end.search(line):
                t, step, task = parse_time(m[1]), int(m[2]), int(m[3])
                send_end[(task, step)] = t

        # Compute initialization duration
        for task in init_start:
            if task in init_end:
                t_init = (init_end[task] - init_start[task]).total_seconds()
                records.append(dict(task_id=task, step=-1, t_init=t_init, t_step=None, t_send=None, worker_id=worker_id))

        # Compute per-step durations
        for (task, step), t0 in step_start.items():
            if (task, step) in step_end:
                t_step = (step_end[(task, step)] - t0).total_seconds()
            else:
                t_step = None
            if (task, step) in send_start and (task, step) in send_end:
                t_send = (send_end[(task, step)] - send_start[(task, step)]).total_seconds()
            else:
                t_send = None
            records.append(dict(task_id=task, step=step, t_init=None, t_step=t_step, t_send=t_send, worker_id=worker_id))

    # --- Build DataFrame ---
    df = pd.DataFrame.from_records(records)

    print(df)

    # --- Aggregate by task ---
    def agg_task(group):
        num_steps = int((group["step"] >= 0).sum())
        t_init = group.loc[group["step"] == -1, "t_init"].sum()
        t_step = group.loc[group["step"] >= 0, "t_step"].sum(skipna=True)
        t_send = group.loc[group["step"] >= 0, "t_send"].sum(skipna=True)
        return pd.Series(dict(num_steps=num_steps, t_init=t_init, t_step=t_step, t_send=t_send, worker_id=group["worker_id"].iloc[0]))

    summary = df.groupby("task_id").apply(agg_task).reset_index().sort_values("task_id")
    # Ensure proper dtypes
    summary = summary.astype({"num_steps": int, "worker_id": int})


    print(summary)
    return summary

def plot_task_times(summary):

    fig, ax = plt.subplots(figsize=(10, 6))

    y_pos = range(len(summary))
    ax.barh(y_pos, summary["t_init"], color='skyblue', label='Initialization')
    ax.barh(y_pos, summary["t_step"], left=summary["t_init"], color='lightgreen', label='Computation')
    ax.barh(y_pos, summary["t_send"], left=summary["t_init"] + summary["t_step"], color='salmon', label='Data Sending')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(summary["task_id"], fontsize=6)
    ax.set_xlabel('Time (seconds)')
    ax.set_title('Task Execution Times')
    ax.legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    summary = parse_logs()
    plot_task_times(summary)