from gurobipy import *

class ExactOAS(object):
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

		self.orders = list(range(totalorders))

		self.setuptime = setuptime


	def modelMIPForOAS(self):
		model = Model('OAS')

		# add variable

		# predecessor
		y = model.addVars([(i,j) for i in self.orders
			for j in self.orders if i != j],
			vtype=GRB.BINARY,
			name='y'
		)

		# accepted order indicator
		x = model.addVars(
			self.orders,
			vtype=GRB.BINARY,
			name='x'
		)

		# completion time
		C = model.addVars(
			self.orders,
			vtype=GRB.CONTINUOUS,
			lb=0,
			name='C'
		)

		# tardiness
		T = model.addVars(
			self.orders, 
			vtype=GRB.CONTINUOUS,
			lb=0,
			name='T'
		)

		# revenue
		R = model.addVars(
			self.orders[1:-1],
			vtype=GRB.CONTINUOUS,
			name='R'
		)

		# add objective
		model.setObjective(sum(R[i] for i in self.orders[1:-1]), sense=GRB.MAXIMIZE)

		# add constraints
		model.addConstrs(
			sum(y[i,j] for j in self.orders[1:] if j != i)\
			== x[i] for i in self.orders[:-1]
		)
		model.addConstrs(
			sum(y[j,i] for j in self.orders[:-1] if j != i)\
			== x[i] for i in self.orders[1:]
		)

		model.addConstrs(
			C[i] + (self.setuptime[i,j]+self.processtime[j])*y[i,j] \
			+ self.deadline[i]*(y[i,j] - 1) <= C[j]
			for i in self.orders[:-1]
			for j in self.orders[1:]
			if i != j
		)

		model.addConstrs(
			(self.releasedate[j] + self.processtime[j])*x[j] \
			+ self.setuptime[i,j]*y[i,j] <= C[j]
			for i in self.orders[:-1]
			for j in self.orders[1:]
			if i != j
		)

		model.addConstrs(
			C[i] <= self.deadline[i]*x[i]
			for i in self.orders
		)

		model.addConstrs(
			T[i] >= C[i] - self.duedate[i]
			for i in self.orders
		)

		model.addConstrs(
			T[i] <= (self.deadline[i] - self.duedate[i])*x[i]
			for i in self.orders
		)

		model.addConstrs(
			R[i] <= self.revenue[i]*x[i] - T[i]*self.weight[i]
			for i in self.orders[1:-1]
		)

		model.addConstrs(
			R[i] >= 0 for i in self.orders[1:-1]
		)

		# for dummy orders
		model.addConstr(C[0] == 0)
		model.addConstr(
			C[self.n+1] == max(
				[self.deadline[i] for i in self.orders[1:-1]]
			)
		)
		model.addConstr(x[0] == 1)
		model.addConstr(x[self.n+1] == 1)

		model.optimize()

		return model