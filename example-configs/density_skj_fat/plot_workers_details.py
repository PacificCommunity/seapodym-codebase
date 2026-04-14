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

    for file in files:
        print(f"Parsing {file}")

        with open(file) as f:

            for line in f:

                # ---------- timestamp ----------
                ts_m = re.search(ts_re, line)
                if not ts_m:
                    continue
                ts = datetime.strptime(ts_m.group(1), "%Y-%m-%d %H:%M:%S.%f")

                # ---------- worker ----------
                w_m = re.search(worker_re, line)
                if not w_m:
                    continue
                worker = int(w_m.group(1))

                # ---------- task ----------
                t_m = re.search(task_re, line)
                if not t_m:
                    continue
                task = int(t_m.group(1))

                # =====================================================
                # 1. CLASSIFY EVENT TYPE FIRST (CRITICAL FIX)
                # =====================================================

                is_init = "initialization of task" in line
                is_put = "send data for step" in line
                is_notify = "notify manager after step" in line
                is_step = ("step" in line and "of task" in line and not is_put and not is_notify and not is_init)

                if is_init:
                    phase = "init"
                    step = None
                    key = (worker, task, "init")

                elif is_step:
                    phase = "step"
                    step = int(re.search(step_re, line).group(1))
                    key = (worker, task, "step", step)

                elif is_put:
                    phase = "put"
                    step = int(re.search(step_re, line).group(1))
                    key = (worker, task, "put", step)

                elif is_notify:
                    phase = "notify"
                    step = int(re.search(step_re, line).group(1))
                    key = (worker, task, "notify", step)

                else:
                    continue

                # ---------- direction ----------
                if ">>>" in line:
                    direction = "start"
                elif "<<<" in line:
                    direction = "end"
                else:
                    continue

                # ---------- match ----------
                if direction == "start":
                    active[key] = ts
                else:
                    if key not in active:
                        continue

                    records.append({
                        "worker_id": worker,
                        "task_id": task,
                        "phase": phase,
                        "step": step,
                        "t_start": active.pop(key),
                        "t_end": ts,
                    })

    df = pd.DataFrame(records)

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

    if not df.empty:
        df = df.sort_values(["worker_id", "t_start"])

    plot_gantt(df)



