'''
'''
from math import ceil
from random import uniform

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

		# revenue load ratios
		self.RLR1 = {}
		self.RLR2 = {}


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
			try:
				nextorderindex = solution.index(i)
			except ValueError:
				print('ValueError occured')
				print(solution)
				print(i, min_, max_)
				raise ValueError()
			nextorder = nextorderindex + 1 # order numbers are 1 to n

			# total time incurred = setup time + process time
			time = self.setuptime[prevorder, nextorder] \
				+ self.processtime[nextorder]

			# absolute time
			completiontimes[nextorderindex] = \
				max(timeelapsed, self.releasedate[nextorder]) + time
			timeelapsed = completiontimes[nextorderindex]

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


	def calculateStartTime(self, solution, completiontimes=None):
		starttimes = [0 for _ in solution]
		if completiontimes == None:
			completiontimes = self.calculateCompletionTimes(solution)

		seq = max(solution)
		while seq > 0:
			ind = solution.index(seq)
			if seq > 1:
				ind_ = solution.index(seq-1)
				setup = self.setuptime[ind_+1, ind+1]
			else:
				setup = self.setuptime[0, ind+1]
			process = self.processtime[ind+1]
			compl = completiontimes[ind]
			starttimes[ind] = compl - process - setup
			seq -= 1
		return starttimes


	def calculatePreemption(self, solution):
		# array to contain postive value if 
		# order processing starts before release date
		preemption = [0 for _ in solution]

		# get start times
		starttimes = self.calculateStartTime(solution)
		for ind, start in enumerate(starttimes, start=1):
			if solution[ind-1] > 0: 
				# order is scheduled
				preemption[ind-1] = self.releasedate[ind] - start
		return preemption


	@staticmethod
	def rejectOrder(solution, indextoreject=None):
		sequencenum = solution[indextoreject]
		newsolution = solution[:]
		newsolution[indextoreject] = 0 #reject order

		# adjust sequence for orders
		newsolution = [
			(lambda x: x-1 if x>sequencenum else x)(x) 
			for x in newsolution]
		return newsolution


	@staticmethod
	def insertOrder(solution, ordertoinsert, insertloc):
		newsolution = solution[:]
		# adjust sequence
		newsolution = [
			(lambda x: x+1 if x >= insertloc else x)(x)
			for x in newsolution
		]
		ordertoinsertind = ordertoinsert - 1
		newsolution[ordertoinsertind] = insertloc
		return newsolution


	def rejectInfeasibleOrders(self, solution, mostnegative=True):
		revenue = self.calculateRevenue(solution)
		while any([revenue[i] < 0 for i, seq in enumerate(solution) if seq > 0]):
			# solution contains scheduled order(s) with nonpositve revenue
			# this will reject orders that might be scheduled but have zero revenue
			if mostnegative:
				# start by rejecting order with min revenue
				indnegativerevenue = revenue.index(min(revenue))
			else:
				# find first non postive
				for ind, rev in enumerate(revenue):
					if rev <=0 and solution[ind] > 0:
						indnegativerevenue = ind
						break
			# reject order
			newsolution = TabuSearchOAS.rejectOrder(solution, indnegativerevenue)

			# print(newsolution)
			revenue = self.calculateRevenue(newsolution)
			solution = newsolution
		# print(solution, sum(revenue))
		# print('preemption: ', self.calculatePreemption(solution))
		return solution, revenue


	def generateRLR1(self):
		for order in range(1, self.n+1):
			avgsetup = sum(self.setuptime[:self.n+1][order])/(self.n + 1)
			self.RLR1[order] = self.revenue[order]/(
				self.processtime[order] + avgsetup)


	def generateRLR2(self, solution):
		self.RLR2 = {}
		prevorder  = 0 # 0 is dummy order
		min_, max_ = max(min(solution), 1), max(solution)

		# iterate over all sequence in solution
		for seq in range(min_, max_+1):
			orderind = solution.index(seq)
			order = orderind + 1
			self.RLR2[order] = self.revenue[order]/(
				self.processtime[order] + \
				self.setuptime[prevorder, order])
			prevorder = order


	def createInitialSolution(self):
		self.generateRLR1()

		sortedorders = {order: ratio 
		for order, ratio in sorted(
			self.RLR1.items(), 
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
		# print(solution, sum(revenue))
		print('Starting with solution: ' +\
			self.displaySolution(solution, sum(revenue))
		)
		
		# iteratively remove infeasible orders
		# causing negative revenue
		solution, revenue = self.rejectInfeasibleOrders(solution)
		# print(solution, sum(revenue))
		print('Initial feasible solution: ' +\
			self.displaySolution(solution, sum(revenue))
		)

		return solution


	def checkFeasibility(self, solution, completiontimes=None):
		feasible = True
		if completiontimes == None:
			# calculate if not available
			completiontimes = self.calculateCompletionTimes(solution)

		starttimes = self.calculateStartTime(solution, completiontimes)

		for ind, seq in enumerate(solution, start=0):
			if seq > 0:
				# scheduled
				# check if completion time is within deadline
				if completiontimes[ind] > self.deadline[ind+1]:
					feasible = False
					break

				# check if start time is past release date
				if starttimes[ind] < self.releasedate[ind+1]:
					feasible = False
					break
			
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
					# same swap occured within tabu tenure
					continue
				elif (order2, order1) in self.tabulist:
					# same swap occured within tabu tenure
					continue
				else:
					# make copy
					neighbor = self.currsol[:]
					# swap
					neighbor[order1-1], neighbor[order2-1] = \
						neighbor[order2-1], neighbor[order1-1]

					# reject infeasible orders
					# neighbor, _ = self.rejectInfeasibleOrders(neighbor)

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


	@staticmethod
	def rouletteWheelSelection(choices):
		max = sum(choices.values())
		pick = uniform(0, max)
		current = 0
		for key, value in choices.items():
			current += value
			if current > pick:
				return key


	def localSearch(self, solution, revenuesum):
		self.generateRLR2(solution)

		# RLR1 for the rejected orders in solution
		rejectedRLR1 = {order: self.RLR1[order]
			for order, seq in enumerate(solution, start=1)
			if seq == 0}

		while len(rejectedRLR1) > 0:
			ordertoinsert = TabuSearchOAS.rouletteWheelSelection(
				rejectedRLR1
			)

			for insertloc in range(1, max(solution) + 2):
				newsolution = TabuSearchOAS.insertOrder(solution, ordertoinsert, insertloc)
				newsolution, newrevenue = self.rejectInfeasibleOrders(newsolution)

				if sum(newrevenue) > 0.998*revenuesum:
					solution = newsolution
					break
				else:
					continue
			del rejectedRLR1[ordertoinsert]
		return solution


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
				print('Improved Solution: ' \
					+ self.displaySolution(self.bestsol, self.bestrev))
				noimprovementcount = 0
			else:
				noimprovementcount += 1
			if noimprovementcount >= 50:
				terminate = True

			# local search procedure for add and insert orders
			for _ in range(ceil(self.n/6)):
				self.currsol = self.localSearch(self.currsol, currsolrev)

		print('Final Solution: ' \
			+ self.displaySolution(self.bestsol, self.bestrev))
		print('Best Solution: ', self.bestsol, self.bestrev)
		completiontimes = self.calculateCompletionTimes(self.bestsol)
		print('Completion Time: ', completiontimes)
		print('Revenue: ', self.calculateRevenue(self.bestsol))
		print('Release: ', self.releasedate[1:-1])
		print('Start  : ', self.calculateStartTime(self.bestsol))
		print('Premp  : ', self.calculatePreemption(self.bestsol))

		return self.bestrev


	def displaySolution(self, solution, revenue=None):
		if revenue == None:
			revenue = sum(self.calculateRevenue(solution))
		return f"Orders Accepted = {sum(s>0 for s in solution)} " + \
			f"Revenue = {revenue}" 