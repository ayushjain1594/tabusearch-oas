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
		completiontime = [0 for _ in solution]

		# calculate completion time for all
		# based on sequence

		# calculate revenue array


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

		# calculate revenue assuming seqeunce 
		# created by sorted orders

		# while revenue_i < 0 for any i
		# 	remove order i
		# 	adjust sequence
		# 	recalculate revenue