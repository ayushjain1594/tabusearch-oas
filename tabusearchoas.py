'''
'''
from math import ceil

class TabuSearchOAS:
	"""
	"""
	def __init__(self, 
		totalorders, 
		releasedate,
		processtime,
		duedate,
		deadline,
		revenue,
		weight,
		setuptime):

		self.n = totalorders - 2 #dummy orders
		self.releasedate = releasedate
		self.processtime = processtime
		self.duedate = duedate
		self.deadline = deadline
		self.revenue = revenue
		self.weight = weight

		self.setuptime = setuptime

		# initial, current, best solution
		self.initsol = None
		self.currsol = None
		self.bestsol = None

		# best revenue
		self.bestrev = 0
		self.tabutenure = ceil(2*self.n/3)

		# tabu list
		self.tabulist = []
		self.tabulistlen = 0


	def calculateCompletionTimes(self, solution):
		# array to hold completion time for each order
		completiontimes = [0 for _ in solution]

		# calculate completion time for all
		# based on sequence
		timeelapsed = 0
		prevorder  = 0 # 0 is dummy order
		min_, max_ = max(min(solution), 1), max(solution)

		# iterate over all sequence in solution
		for i in range(min_, max_+1):
			nextorderindex = solution.index(i)
			nextorder = nextorderindex + 1 # order numbers are 1 to n

			# total time incurred = setup time + process time
			time = self.setuptime[prevorder, nextorder] \
				+ self.processtime[nextorder]

			# absolute time
			completiontimes[nextorderindex] = timeelapsed + time
			timeelapsed = timeelapsed + time

			prevorder = nextorder

		# print(completiontimes, max(completiontimes))
		return completiontimes


	def calculateRevenue(self, solution):
		# get completion time for the solution
		completiontimes = self.calculateCompletionTimes(solution)

		# measure of delay in completing order
		# enumerate starting 1
		tardiness = [max(0, ctime - self.duedate[ind])
			for ind, ctime in enumerate(completiontimes, start=1)]
		# print(tardiness)

		revenue = [
			self.revenue[ind]*(1 if solution[ind-1] > 0 else 0) \
			- self.weight[ind]*delay
			for ind, delay in enumerate(tardiness, start=1)]

		# print(revenue, sum(revenue), "\n")
		return revenue


	def rejectInfeasibleOrders(self, solution, mostnegative=True):
		revenue = self.calculateRevenue(solution)
		while any([revenue[i] <= 0 for i, seq in enumerate(solution) if seq > 0]):
			# solution contains scheduled order(s) with nonpositve revenue
			# this will reject orders that might be scheduled but have zero revenue
			if mostnegative:
				# start by rejecting order with min revenue
				indmaxnegativerevenue = revenue.index(min(revenue))
			else:
				# find first non postive
				for ind, rev in enumerate(revenue):
					if rev <=0 and solution[ind] > 0:
						indmaxnegativerevenue = ind
						break
			# print(indmaxnegativerevenue)
			sequencenum = solution[indmaxnegativerevenue]

			# copy solution
			newsolution = solution[:]
			newsolution[indmaxnegativerevenue] = 0 # reject order

			# adjust sequence after rejecting order
			newsolution = [
				(lambda x: x-1 if x>sequencenum else x)(x) 
				for x in newsolution]

			# print(newsolution)
			revenue = self.calculateRevenue(newsolution)
			solution = newsolution
			# print(solution, sum(revenue))
		return solution, revenue


	def createInitialSolution(self):
		revenueloadratio = {}

		for i in range(1, self.n+1):
			avgsetup = sum(self.setuptime[:][i-1])/(self.n + 1)
			revenueloadratio[i] = self.revenue[i]/(
				self.processtime[i] + avgsetup)

		sortedorders = {order: ratio 
		for order, ratio in sorted(
			revenueloadratio.items(), 
			key=lambda item: item[1],
			reverse=True
			)
		}

		# array containing seq number of order 1 to n
		# ith index -> seq of order i+1
		# 0 if order is rejected
		solution = [0]*self.n
		seq = 1
		# solution with seq based on ratio
		for order, _ in sortedorders.items():
			if self.revenue[order] <= 0:
				# for orders with no revenue regardless
				# of completing within time
				continue
			solution[order-1] = seq
			seq += 1

		revenue = self.calculateRevenue(solution)
		print(solution, sum(revenue))
		
		# iteratively remove infeasible orders
		# causing negative revenue
		solution, revenue = self.rejectInfeasibleOrders(solution)
		print(solution, sum(revenue))

		return solution


	def checkFeasibility(self, solution, completiontimes=None):
		feasible = True
		if completiontimes == None:
			# calculate if not available
			completiontimes = self.calculateCompletionTimes(solution)

		for orderno, ctime in enumerate(completiontimes, start=1):
			# check if completion time is within due date
			if ctime > self.duedate[orderno]:
				feasible = False
				break

			# check if release date could cause infeasibility
			# INSERT LOGIC HERE
		return feasible


	def generateNeighbors(self):
		# swap & insertion
		# insertion here refers to including one order 
		# and rejecting another when one of the seq swapped is 0
		neighbors = []
		swaps = {}
		for order1, seq1 in enumerate(self.currsol[:-1], start=1):
			for order2, seq2 in enumerate(self.currsol[order1:], start=order1+1):
				if seq1 == 0 and seq2 == 0:
					# same neighbor
					continue
				elif (order1, order2) in self.tabulist:
					continue
				elif (order2, order1) in self.tabulist:
					continue
				else:
					# make copy
					neighbor = self.currsol[:]
					# swap
					neighbor[order1-1], neighbor[order2-1] = \
						neighbor[order2-1], neighbor[order1-1]

					# reject infeasible orders
					neighbor, _ = self.rejectInfeasibleOrders(neighbor)

					# track swap
					swaps[tuple(neighbor)] = [(order1, order2), (order2, order1)]
					# append to neighbors
					neighbors.append(neighbor)
		return neighbors, swaps


	def findBestNeighbor(self, neighbors):
		bestneighborrev = -99999
		bestneighbor = None
		for neighbor in neighbors:
			revenue = sum(self.calculateRevenue(neighbor))
			if revenue > bestneighborrev:
				bestneighborrev = revenue
				bestneighbor = neighbor

		return bestneighbor


	def updateTabuList(self, swaps):
		if self.tabulistlen >= self.tabutenure:
			# remove first two moves
			self.tabulist = self.tabulist[2:]
			self.tabulist.extend(swaps)
		else:
			self.tabulist.extend(swaps)
			self.tabulistlen += len(swaps)


	def tabuSearchAlgorithm(self, threshould=50):

		# create initial solution
		self.initsol = self.createInitialSolution()
		self.currsol = self.initsol[:]

		terminate = False
		# count of iterations with no improvement in sol
		noimprovementcount = 0
		print("Starting Tabu Search")
		# iterate until termination criteria is met
		while not terminate:
			neighbors, swaps = self.generateNeighbors()

			# best solution among neighbors
			bestneighbor = self.findBestNeighbor(neighbors)
			self.currsol = bestneighbor

			# update the tabu list
			self.updateTabuList(swaps.get(tuple(bestneighbor), []))

			currsolrev = sum(self.calculateRevenue(self.currsol))

			if currsolrev > self.bestrev:
				self.bestsol = self.currsol
				self.bestrev = currsolrev
				print(f'Improved Solution {self.bestsol} and Reveue {self.bestrev}')
				noimprovementcount = 0
			else:
				noimprovementcount += 1
			if noimprovementcount >= 50:
				terminate = True

			# write here randomized local search procedure

		print(noimprovementcount, self.bestsol, self.bestrev)
		return self.bestsol