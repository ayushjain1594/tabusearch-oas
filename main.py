'''
'''

from numpy import zeros
from pathlib import Path

from tabusearchoas import TabuSearchOAS

def fetchDataFromLines(line):
	return [float(data) for data in line.replace('\n', '').split(',')]


def fetchData(filepath):
	data = open(Path.cwd().joinpath(filepath), 'r')
	lines = data.readlines()

	releasedate = fetchDataFromLines(lines[0])
	processtime = fetchDataFromLines(lines[1])
	duedate = fetchDataFromLines(lines[2])
	deadline = fetchDataFromLines(lines[3])
	revenue = fetchDataFromLines(lines[4])
	weight = fetchDataFromLines(lines[5])

	ordercount = len(releasedate)

	setuptime = zeros([ordercount, ordercount], dtype=float)
	for ind in range(ordercount):
		setuptime[ind] = fetchDataFromLines(lines[ind+6])

	# print(releasedate)
	# print(processtime)
	# print(duedate)
	# print(deadline)
	# print(revenue)
	# print(weight)

	# print(setuptime)
	return (ordercount, releasedate, processtime, 
		duedate, deadline, revenue,
		weight, setuptime)


def solveOAS(args):
	ts = TabuSearchOAS(*args)
	# ts.createInitialSolution()
	ts.tabuSearchAlgorithm()

if __name__ == '__main__':
	args = fetchData("Dataset_OAS/Dataset_OAS/100orders/Tao9/R9/Dataslack_100orders_Tao9R9_9.txt")
	solveOAS(args)