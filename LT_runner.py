from simulation import Well, LTSimulation, Cell
import pickle
import pprint
import time, math
'''
EventLoop is the bread and butter of this simulation, lets take a look at it.
It's initilized as an empty dict in LTSim, and each cell is given a reference to it
The format is :
	Events = {1122: <--tick number
				{
					cellid_1 : 'division',
					cellid_2 : 'edit',
					cellid_3 : 'division',
					...
				}
				...
			}
'''

# Constants
SPLIT_TIMES = [1440, 2880]# 4320]
DUMP_FILE = f"./simdump_{math.floor(time.time())}"

def main():
	s = LTSimulation({'split_times':SPLIT_TIMES})
	print("\n***********************************************************************")
	print("** Starting Lineage Tracing Simuation with the following parameters: **\n")
	pprint.PrettyPrinter().pprint(s.parameters)
	print("")
	print("***********************************************************************\n")

	print("***********************************************************************")
	print("********************* Pre-Simulation Expansion ************************")
	print("***********************************************************************\n")

	while s.getActiveCells() < s.parameters["initialCellCount"]*s.parameters["concurrentLineages"]:
		try:
			s.init[s.currentTick]
		except KeyError:
			s.currentTick += 1
			continue
		else:
			for cell in s.init[s.currentTick]:
				s.divideCell(cell, s.currentTick, init=True)
			s.currentTick += 1

	for tick in [t for t in s.init if int(t) >= s.currentTick]:
		dTime = int(tick) - s.currentTick
		for c in s.init[tick]:
			cell = s.cells[c]
			if not cell.removed:
				cell.createEvent(dTime, "division")

	print("***********************************************************************")
	print(f"*************** Simulation Beginning with {s.getActiveCells()} cells *****************")
	print("***********************************************************************")
	
	s.populateRoots()
	s.currentTick = 0
	for tick in range(s.getMaxTicks()): # Main Simulation Loop
		s.currentTick = tick
		try:
			s.events[tick] # Are there any events at this tick?
		except KeyError:
			continue #catch the error and move on.
		else:	
			for cell in s.events[tick]:
				if not cell == 'SPLIT':
					if not s.cells[cell].removed:
						event = s.events[tick][cell]
						if event == 'division':
							# print(f"Dividing cell: {cell}")
							s.divideCell(cell, tick, init=False)
						else:
							raise Exception("You fucked up.")
				else:
					print(f"Splitting {s.getActiveCells()} cells.")
					s.splitCells(s.events[tick]['SPLIT'])

	s.exportData(DUMP_FILE);
	total_len = 0
	total_ins = 0
	have = 0
	have_not = 0
	for c in s.cells:
		cell = s.cells[c]
		if not cell.removed:
			if cell.ins:
				total_len += len("".join(cell.ins))
				total_ins += 1
				have += 1
			else:
				have_not += 1

	print(f"Average insertion length = {total_len/total_ins}")
	print(f"Number of active cells at end of simulation: {s.getActiveCells()}")
	print(f"Cells with insertions: {have}")
	print(f"Cells without insertions: {have_not}")	

if __name__ == '__main__':
	main()