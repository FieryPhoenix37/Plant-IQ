"""Quick demo: build the primary-school network and run the authors' SEIRS+ K-12 sim loop.

Parameters here are placeholders to show the code runs; they are not the paper's calibrated values.
Run from the repo root:  python scripts/demo_k12_sim.py
"""

import warnings

import numpy as np

from school.seirsplus_k12.models import ExtSEIRSNetworkModel
from school.seirsplus_k12.networks import generate_K5_school_contact_network
from school.seirsplus_k12.sim_loops import run_tti_sim

np.random.seed(1)
nets, info = generate_K5_school_contact_network(
    num_grades=6, num_classrooms_per_grade=4, class_sizes=20, num_student_blocks=2, num_staff=24
)
G = nets["onsite-all"]
print(f"network: {G.number_of_nodes()} nodes "
      f"({len(info['studentIDs'])} students, {len(info['teacherIDs'])} teachers, {len(info['staffIDs'])} staff)")

model = ExtSEIRSNetworkModel(G=G, p=0.2, beta=0.3, sigma=1 / 3, lamda=1 / 2.2, gamma=1 / 4, initE=1,
                             transition_mode="time_in_state", isolation_time=10)
week = [nets["onsite-all"]] * 5 + [nets["offsite-all"]] * 2  # school Mon-Fri, no school contacts on weekends

warnings.simplefilter("ignore")
run_tti_sim(model, T=60, max_dt=0.5, average_introductions_per_day=1 / 7, parameter_cadences={"G": week})

ever_infected = int(model.numNodes - model.numS[model.tidx])
print(f"after {model.t:.0f} days: {ever_infected} people ever infected")
