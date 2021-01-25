'''
'''

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


	def calculateRevenue(self, solution):

		# array to hold completion time for each order
		completiontime = [0 for _ in solution]

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
			completiontime[nextorderindex] = timeelapsed + time
			timeelapsed = timeelapsed + time

			prevorder = nextorder

		# print(completiontime, max(completiontime))

		# measure of delay in completing order
		# enumerate starting 1
		tardiness = [max(0, ctime - self.duedate[ind])
			for ind, ctime in enumerate(completiontime, start=1)]
		print(tardiness)

		revenue = [
			self.revenue[ind]*(1 if solution[ind-1] > 0 else 0) \
			- self.weight[ind]*delay
			for ind, delay in enumerate(tardiness, start=1)]

		print(revenue, sum(revenue), "\n")
		return revenue


	def createInitialSolution(self, mostnegative=True):
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

		print(solution)
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
			print(indmaxnegativerevenue)
			sequencenum = solution[indmaxnegativerevenue]

			# copy solution
			newsolution = solution
			newsolution[indmaxnegativerevenue] = 0 # reject order

			# adjust sequence after rejecting order
			newsolution = [
				(lambda x: x-1 if x>sequencenum else x)(x) 
				for x in newsolution]

			print(newsolution)
			revenue = self.calculateRevenue(newsolution)
			solution = newsolution

		return solution