import re
import glob
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

# ---------- CONFIG ----------
LOG_PATTERN = "log_taskfunc*.txt"

COLORS = {
    "init": "green",
    "step": "lightblue",
    "put": "red",
    "notify": "orange",
}

# ---------- PARSER ----------
def parse_logs(file_pattern):
    records = []
    active = {}

    files = sorted(glob.glob(file_pattern))
    print(f"Found {len(files)} log files")

    timestamp_re = r"\[(.*?)\]"
    worker_re = r"\[log(\d+)\]"
    task_re = r"task id (\d+)"
    step_re = r"step (\d+)"

    for file_path in files:
        print(f"Parsing {file_path}")

        with open(file_path) as f:
            for line in f:
                # --- Timestamp ---
                ts_match = re.search(timestamp_re, line)
                if not ts_match:
                    continue
                timestamp = datetime.strptime(
                    ts_match.group(1),
                    "%Y-%m-%d %H:%M:%S.%f"
                )

                # --- Worker ID ---
                worker_match = re.search(worker_re, line)
                if not worker_match:
                    continue
                worker_id = int(worker_match.group(1))

                # --- Direction ---
                if ">>>" in line:
                    direction = "start"
                elif "<<<" in line:
                    direction = "end"
                else:
                    continue

                # --- Task ID ---
                task_match = re.search(task_re, line)
                if not task_match:
                    continue
                task_id = int(task_match.group(1))

                # --- Phase detection ---
                if "initialization" in line:
                    phase = "init"
                    step = None
                elif "send data" in line:
                    phase = "put"
                    step = int(re.search(step_re, line).group(1))
                elif "notify manager" in line:
                    phase = "notify"
                    step = int(re.search(step_re, line).group(1))
                elif ">>> step" in line or "<<< step" in line:
                    phase = "step"
                    step = int(re.search(step_re, line).group(1))
                else:
                    continue

                key = (worker_id, task_id, phase, step)

                if direction == "start":
                    active[key] = timestamp
                else:
                    if key in active:
                        t_start = active.pop(key)
                        t_end = timestamp

                        records.append({
                            "worker_id": worker_id,
                            "task_id": task_id,
                            "phase": phase,
                            "step": step,
                            "t_start": t_start,
                            "t_end": t_end,
                            "source_file": file_path,  # 👈 useful debug
                        })

    df = pd.DataFrame(records)

    if df.empty:
        print("Warning: No records parsed!")
    else:
        print(f"Parsed {len(df)} intervals")

    return df


# ---------- GANTT PLOT ----------
def plot_gantt(df):
    fig, ax = plt.subplots(figsize=(14, 6))

    # Normalize time
    t0 = df["t_start"].min()
    df["start_s"] = (df["t_start"] - t0).dt.total_seconds()
    df["end_s"] = (df["t_end"] - t0).dt.total_seconds()

    # Plot
    for _, row in df.iterrows():
        ax.barh(
            y=row["worker_id"],
            width=row["end_s"] - row["start_s"],
            left=row["start_s"],
            color=COLORS.get(row["phase"], "gray"),
        )

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Worker ID")
    ax.set_title("MPI Worker Timeline (All Logs)")

    # Clean y-axis
    workers = sorted(df["worker_id"].unique())
    ax.set_yticks(workers)

    # Legend
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=c)
        for c in COLORS.values()
    ]
    labels = list(COLORS.keys())
    ax.legend(handles, labels)

    plt.tight_layout()
    plt.show()


# ---------- MAIN ----------
if __name__ == "__main__":
    df = parse_logs(LOG_PATTERN)

    # Sort for nicer plotting
    df = df.sort_values(["worker_id", "t_start"])

    print(df.head())

    plot_gantt(df)

