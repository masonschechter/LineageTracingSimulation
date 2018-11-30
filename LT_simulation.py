import math
import numpy as np
import pickle
import pprint
import math
from itertools import product

INSERTION_LEN_RANGE = [0, 1, 2, 3, 4, 5, 6, 7]
INSERTION_LEN_FREQ = [.05, .15, .15, .25, .15, .1, .1, .05]
NUCLEOTIDES = ["G", "C", "A", "T"]
NUCLEOTIDE_BIAS = [.35, .40, .15, .10]

#***********************************#
#*********** CELL CLASS ************#
#***********************************#
class Cell():
	
	def __init__(self, well = "", parent = None, ins = [], createdAt = 0, eventLoop = {}, initLoop = {}, init=False):
		self.id = "".join(list(np.random.choice(["G", "C", "A", "T"], size=20))) #TODO Check uniqueness
		self.well = well
		self.ins = ins
		self.parent = parent
		self.daughters = []
		self.createdAt = createdAt
		self.events = eventLoop #main simulation event loop
		self.initLoop = initLoop #pre-simulation population growth loop
		self.dividedAt = -1
		self.removed = False
		self.init = init
		if not self.init:
			self.future = self.predictFuture()
		else:
			self.future = self.predictFuture(init=True)
	
	# TODO  	
	def predictFuture(self, init=False):
		dTime = self.createdAt + int(np.random.normal(24*60,120)) #normal distribution of HEK293T doubling time, mean: 24 hrs		
		if not init:
			self.createEvent(dTime,"division")
		else: 
			self.createEvent(dTime,"division", init=True)


	def divide(self, divisionTime, init=False):
		# Division occurs in two steps:
		# 1) Create two new cell objects as the daughter cells, initialize them with the current dtime as their starting time
		# 2) append those cells to this cells daughters array and return them so the LT sim can also keep track of them
		if not init:
			d1 = Cell(self.well, self, self.ins[:], divisionTime+1, self.events, self.initLoop)
			d2 = Cell(self.well, self, self.ins[:], divisionTime+1, self.events, self.initLoop)
		else:
			d1 = Cell(self.well, self, self.ins[:], divisionTime+1, self.events, self.initLoop, init=True)
			d2 = Cell(self.well, self, self.ins[:], divisionTime+1, self.events, self.initLoop, init=True)
		self.daughters.append(d1)
		self.daughters.append(d2)
		return d1,d2

	def addInsertion(self, seq):
		self.ins.append(seq)


	def createEvent(self, time, event, init=False):
		if not init: # we're in the main loop of the simulation
			try:
				self.events[time][self.id] = event #there's already an event at this time
			except:
				self.events[time] = {self.id:event} #if not, create it.
		else: # we're in the pre-simulation population expansion
			try:
				self.initLoop[time][self.id] = event
			except:
				self.initLoop[time] = {self.id:event}

	# def addWell(self, wellNum):
		# self.well += wellNum

	def __str__(self):
		return self.id

#***********************************#
#******** SIMULATION CLASS *********#
#***********************************#
class LTSimulation():

	def __init__(self, parameters=None):
		self.parameters = {
							"tickTime":1, #1 minute Ticks
							"splits":2,
							"initialCellCount": 10000,
							"totalTime": 0,
							"totalTicks":0,
							"concurrentLineages":4,
							"splitSize":3
							}
		if parameters:
			for parameter in parameters:
				self.parameters[parameter] = parameters[parameter]
		self.parameters["totalTime"] = int(60*24*(self.parameters["splits"]+2))
		self.parameters["totalTicks"] = int(self.parameters["totalTime"]/self.parameters["tickTime"])
		self.removedCells = {}
		self.events = {}
		self.init = {}
		self.initialCell = Cell("init", None, [], 0, self.events, self.init, init=True)
		self.cells = {self.initialCell.id:self.initialCell}
		# self.rootWells = [Well(str(x), None) for x in range(self.concurrentLineages)]
		for split in self.parameters['split_times']:
			self.events[split] = {'SPLIT':self.parameters['split_times'].index(split)+1}
		# self.wells = {"0":self.rootWell}
		# self.createWells()
		self.currentTick = 0

	# def createWells(self):
	# 	q = []
	# 	q.append(x for x in self.rootWells)
	# 	while q:
	# 		well = q.pop(0)
	# 		d1 = Well(well.id+"0", well)
	# 		d2 = Well(well.id+"1", well)
	# 		well.children[d1.id] = d1
	# 		well.children[d2.id] = d2
	# 		self.wells[d1.id] = d1
	# 		self.wells[d2.id] = d2
	# 		if len(well.id) < self.parameters["splits"]:
	# 			q.append(d1)
	# 			q.append(d2)

	def divideCell(self, cellid, divisionTime, init=False): # **only function with cell.id as argument, all others have cell object**
		cell = self.cells[cellid]
		cell.removed = True
		cell.dividedAt = divisionTime
		if not init:
			d1,d2 = cell.divide(divisionTime)
			self.editCell(d1)
			self.editCell(d2)
		else:
			d1,d2 = cell.divide(divisionTime, init=True)
		self.cells[d1.id] = d1
		self.cells[d2.id] = d2
		# self.wells[d1.well].cells[d1.id] = d1
		# self.wells[d2.well].cells[d2.id] = d2

	def editCell(self,cell):
		'''	1) Determine if there's an edit
				1a)	If there's an insertion, (2)
				1b) If there's a deletion, mark the cell as 'removed'
			2) Determine the insertion length based on the ins len frequency
			3) Determine the nucleotides that get added based on TdT's bais.'''

		insLength = len("".join(cell.ins))
		eChance = 50 if not insLength else 75/(.5*insLength) #chance of a cut as a function of insertion length (stgRNA length)
		eDist = np.random.normal(eChance,math.sqrt(eChance))/100 #normal distribution of edit chance, mean of eChance, to determine if an edit occurs
		
		if eDist <= 0: ## adjusting for edge cases which would make the downstream insertion probability negative
			eDist = 0
		elif eDist >= 1:
			eDist = 1

		if np.random.choice((0,1),p=[1-eDist,eDist]): #there's an edit
			if np.random.choice((0,1), p=[0.2,0.8]): #there's an insert
				seq = ''
				for x in range(np.random.choice(INSERTION_LEN_RANGE, p=INSERTION_LEN_FREQ)): #determine the insertion length
					seq += np.random.choice(NUCLEOTIDES, p=NUCLEOTIDE_BIAS) #insertion sequence
				cell.addInsertion(seq)
			else: #Deletion
				cell.removed = True

	def populateRoots(self): #split initial cells into each lineage's root well
		for cell in self.cells:
			cell = self.cells[cell]
			well = np.random.choice(range(self.parameters["concurrentLineages"]))
			cell.well = str(well)

	def splitCells(self, splitNumber):
		for cellId in self.cells:
			cell = self.cells[cellId]
			if not cell.removed:
				cell.well += str(np.random.choice(range(self.parameters["splitSize"])))

	def exportData(self, fileName):
		wellMapping = {}
		i = 1
		for x in range(self.parameters["concurrentLineages"]):
			for y, z in product(range(self.parameters["splitSize"]), range(self.parameters["splitSize"])):
				wellMapping[str(x)+str(y)+str(z)] = f'well_{x}{y}{z}'
				i += 1
		insDict = {}
		for c in self.cells:
			cell = self.cells[c]
			if not cell.removed and cell.ins:
				ins = "".join(cell.ins)
				if cell.well in wellMapping:
					well = wellMapping[cell.well]
					if ins not in insDict:
						insDict[ins] = {"umis":{well:[]}, "counts":{well:1}}
						insDict[ins]["umis"][well].append(cell.id)
					else:
						if well not in insDict[ins]["umis"] and well not in insDict[ins]["counts"]:
							insDict[ins]["umis"][well] = []
							insDict[ins]["umis"][well].append(cell.id)
							insDict[ins]["counts"][well] = 1
						else:
							insDict[ins]["umis"][well].append(cell.id)
							insDict[ins]["counts"][well] += 1
		with open(fileName,'wb') as fname:
			print("Saving to file....")
			pickle.dump(insDict ,fname)

	def printCells(self):
		for cell in self.cells:
			print(cell)

	def getSimulationTime(self):
		return self.parameters["totalTime"]

	def getTickTime(self):
		self.parameters["tickTime"]

	def getActiveCells(self):
		return len([x for x in self.cells if not self.cells[x].removed])

	def getTotalCells(self):
		return len(self.cells)

	def getMaxTicks(self):
		return self.parameters["totalTicks"]

#***********************************#
#******** 	  WELL CLASS   *********#
#***********************************#
class Well():
	def __init__(self, name, parent):
		self.id = name
		self.parent = parent
		self.children = {}
		self.cells = {}


