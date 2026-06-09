import re
import glob
import defopt
import csv


def main(*, dir : str='results', output : str='results/timings.csv'):
    """
    Assemble timings into a csv file
    :param dir directory
    """
    pattern = re.compile(r"time manager\s*=\s*([0-9.]+)")
    rows = []

    for f in glob.glob(dir + '/run*.out'):
        
         print(f'working on {f}...')
        
         parts = f.split("/")[-1].replace(".out","").split("_")
         n = int(parts[1])
         r = int(parts[2])
         host_branch = parts[3]
         host, branch = host_branch.split('-')

         jobid = None
         manager_time = None

         with open(f) as fh:
            for line in fh:
                if line.startswith("jobid="):
                    jobid = line.split("=")[1].strip()
                m = pattern.search(line)
                if m:
                    manager_time = float(m.group(1))
                    print(f'manager time = {manager_time}')

            rows.append((n, jobid, host, branch, manager_time))

         rows.sort()

         with open(output, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["nprocs","jobid","host","branch", "manager_time"])
            writer.writerows(rows)
            
if __name__ == '__main__':
	defopt.run(main)

