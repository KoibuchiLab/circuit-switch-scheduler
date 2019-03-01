# circuit-switch-scheduler
This repo contains the work on job scheduling over circuit-switched network.
## Source Files
### Makefile
This file produces one executable file:
* css.out

    circuit switch scheduler (see details in circuit-switch-scheduler.cc)

Usage: 
> cat workload.txt | ./css.out -T [0-5] -a $node_num (other parameters and usages are the same as cst.out)

### circuit-switch-scheduler.cc
This program schedules the jobs (including flows and pairs) described in the workload file.

#### Example of workload file (e.g., workload.txt)
submit_time run_time node_num source destination flow_id job_id

    1 1 4 1 3 0 0
    1 1 4 2 4 1 0
    1 1 4 2 1 1 0
    1 1 4 3 4 2 0
    5 6 2 0 3 0 1
    7 3 8 5 6 0 2
    7 3 8 4 1 1 2
    7 3 8 3 5 2 2
    7 3 8 2 7 2 2
    7 3 8 6 8 2 2
    9 5 4 8 1 0 3
    9 5 4 2 9 0 3

* Simulation results are saved in output/, sorted by scheduling clocks (e.g., t1, t5, t7, t9,...).
* For each scheduling clock (showed in folder output/t*), there is one system log file (stat) and routing table for each switch (sw*).