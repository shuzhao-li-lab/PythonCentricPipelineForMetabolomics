import csv
import multiprocessing as mp
import sys
import subprocess

def process_job(job):
    subprocess.run(["/Users/mitchjo/Projects/Atlas_10_12_23/process.sh", job["Sequence"], job["Moniker"], job["Filter"], job["MS2_Path"]])

def main():
    jobs = []
    with open(sys.argv[1]) as jobs_csv:
        for entry in csv.DictReader(jobs_csv):
            jobs.append(entry)


    pool = mp.Pool(4)
    pool.map(process_job, jobs)

if __name__ == '__main__':
    mp.freeze_support()
    main()