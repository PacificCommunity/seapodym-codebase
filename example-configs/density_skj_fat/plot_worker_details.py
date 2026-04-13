import pandas as pd
import re
from datetime import datetime
from glob import glob
import os
from pathlib import Path
import defopt


# --- Regex patterns ---
re_task_start = re.compile(
        r"\[(.*?)\].*> task id (\d+) for steps (\d+) to (\d+)"
    )
re_task_end = re.compile(
        r"\[(.*?)\].*< task id (\d+) for steps (\d+) to (\d+)"
    )

# --- Helper to parse timestamps ---
def parse_time(s):
    return datetime.strptime(s, "%Y-%m-%d %H:%M:%S.%f")


def parse_logs(dir):


    # --- Collect raw records ---
    records = []
    all_times = []

    for fname in sorted(dir.glob("log_taskfunc*.txt")):

        worker_match = re.search(r"log_taskfunc(\d+)", os.path.basename(fname))
        worker_id = int(worker_match.group(1)) if worker_match else None

        with open(fname) as f:
            lines = f.readlines()

        task_start, task_end = {}, {}

        for line in lines:
            if m := re_task_start.search(line):
                t, task_id, step_beg, step_end = (
                    parse_time(m[1]),
                    int(m[2]),
                    int(m[3]),
                    int(m[4]),
                )
                task_start[task_id] = dict(
                    t_start=t,
                    step_beg=step_beg,
                    step_end=step_end,
                )
                all_times.append(t)
            elif m := re_task_end.search(line):
                t, task_id, step_beg, step_end = (
                    parse_time(m[1]),
                    int(m[2]),
                    int(m[3]),
                    int(m[4]),
                )
                task_end[task_id] = dict(
                    t_end=t,
                    step_beg=step_beg,
                    step_end=step_end,
                )
                all_times.append(t)

        for task_id in task_start:
            if task_id in task_end:
                beg = task_start[task_id]["step_beg"]
                end = task_start[task_id]["step_end"]
                num_steps = end - beg
                t_start = task_start[task_id]["t_start"]
                t_end = task_end[task_id]["t_end"]
                records.append(
                    dict(
                        worker_id=worker_id,
                        task_id=task_id,
                        num_steps=num_steps,
                        t_start=t_start,
                        t_end=t_end,
                    )
                )

    # --- Convert to DataFrame ---
    df = pd.DataFrame.from_records(records).sort_values(["worker_id", "task_id"]).reset_index(drop=True)

    # --- Convert times to seconds since first timestamp ---
    if not df.empty:
        t0 = min(all_times)
        df["t_start"] = (df["t_start"] - t0).dt.total_seconds()
        df["t_end"] = (df["t_end"] - t0).dt.total_seconds()

    print(df)
    return df

def plot_task_times(df):

    import matplotlib.pyplot as plt

    # Compute task durations
    df["duration"] = df["t_end"] - df["t_start"]

    plt.figure(figsize=(12, 6))

    # Plot a horizontal bar for each task
    for idx, row in df.iterrows():
        plt.barh(
            y=row["worker_id"],            # vertical position = worker
            width=row["duration"],         # horizontal width = task duration
            left=row["t_start"],           # horizontal start = t_start
            height=0.4,                    # bar thickness
            align="center",
            color="skyblue",
            edgecolor="black"
        )
        # Optional: annotate task_id on the bar
        plt.text(
            x=row["t_start"] + row["duration"]/2,
            y=row["worker_id"],
            s=f'{int(row["task_id"])}', # str(row["task_id"]),
            ha="center",
            va="center",
            fontsize=8,
            color="black"
        )

    plt.xlabel("Time (seconds since first task)")
    plt.ylabel("Worker ID")
    plt.title("Task Execution Timeline per Worker")
    plt.yticks(sorted(df["worker_id"].unique()))
    plt.grid(axis="x", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.show()
    
def main(*, dir: Path='.'):
    """
    dir: directory containing the log files
    """
    df = parse_logs(dir)
    plot_task_times(df)


if __name__ == "__main__":
    defopt.run(main)
