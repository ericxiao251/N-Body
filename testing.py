import numpy as np
import pandas as pd
from math import pow, sqrt

def get_position(i):
	return (positions_x[i] + h * velocities_x[i], positions_y[i] + h * velocities_y[i])

def get_velocity(i):
	return (velocities_x[i] + h * (forces_x[i] / masses[i]), velocities_y[i] + h * (forces_y[i] / masses[i]))

def get_forces():
	for q in range(number_of_particles):
		for k in range(number_of_particles):
			if q == k: continue

			x1, y1 = get_position(q)
			x2, y2 = get_position(k)
			dx = x1 - x2
			dy = y1 - y2

			m1, m2 = masses[q], masses[k]
			r = sqrt(pow(dx, 2) + pow(dy, 2))

			if r < EPISLON:
				f = -G * m1 * m2
				F = f / pow(EPISLON, 3)
			else:
				f = -G * (m1 * m2)
				F = f / pow(r, 3)

			forces_x[q] += F * (dx)
			forces_y[q] += F * (dy)

def compute():
	get_forces()
	for i in range(number_of_particles):
		velocities_x[i], velocities_y[i] = get_velocity(i)
		positions_x[i], positions_y[i] = get_position(i)

def run():
	with open('data/theoretical.csv', 'w') as f:
		header = "t,p_id,mass,p_x,p_y,v_x,v_y\n"
		# f.write(header)

		for t in range(time + 1):
			for p in range(number_of_particles):
				results = "{:.6f},{:d},{:.6f},{:.6f},{:.6f},{:.6f},{:.6f}\n"\
				.format(t * h, p, masses[p],
					positions_x[p], positions_y[p],
					velocities_x[p], velocities_y[p]
				)
				f.write(results)
				# print(s)

			compute()

def compare():
	difference = []
	with open("data/theoretical.csv", "r") as file1, open("data/actual.csv", "r") as file2:
		for line1, line2 in zip(file1, file2):
			if line1 != line2:
				difference.append("actual:      {}".format(line1))
				difference.append("theoretical: {}".format(line2))

	with open('data/difference.csv', 'w') as f:
	    for line in difference:
	        f.write(line)

if __name__ == '__main__':
	G = 6.673 * pow(10, -11)
	EPISLON = 0.000000000000000222

	# change these up values
	number_of_particles = 30
	time = 5
	h = 0.05

	header = ["t","p_id","mass","p_x","p_y","v_x","v_y"]
	df = pd.read_csv('data/actual.csv', header=None, names=header)

	masses = df['mass'].tolist()[:number_of_particles]
	velocities_x = df['v_x'].tolist()[:number_of_particles]
	velocities_y = df['v_y'].tolist()[:number_of_particles]
	positions_x = df['p_x'].tolist()[:number_of_particles]
	positions_y = df['p_y'].tolist()[:number_of_particles]
	forces_x = [0] * number_of_particles
	forces_y = [0] * number_of_particles

	run()
	compare()
