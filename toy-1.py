import json
import dimod
import neal
from dwave.system import DWaveSampler, EmbeddingComposite

with open('toy-1-data.json') as f:
    data = json.load(f)

time_to_remove_barrel = data["warehouse"]["timings"]["remove"]
row_capacity = data["warehouse"]["rowCapacity"]
rows = data["warehouse"]["rows"]

# do the picking

picklist = data["picklist"]
removed_barrels = []

for row in rows:
    for b in range(len(row)-1, -1, -1):
        barrel = row[b]
        if barrel["label"] not in picklist: continue
        for rb in range(b+1, len(row)):
            removed_barrels.append(row.pop())
        row.pop()
            
print("data", data)
print("removed", removed_barrels)

# optimise the return of removed barrels
#allocations = [dimod.Integer(label = f"{barrel['label']}",
#                            lower_bound = 0,
#                            upper_bound = len(rows)) 
#                            for barrel in removed_barrels]

def removal_time_cost(barrel, row):
    barrel_next_use = barrel["nextUse"]
    time_cost = 0
    for b in range(len(row)):
        row_barrel_next_use = row[b]["nextUse"]
        if row_barrel_next_use <= barrel_next_use:
            if row_barrel_next_use < barrel_next_use:    
                time_cost += time_to_remove_barrel # * (len(row) - b)
            # if equal then no additional cost, row needs moving anyway
            break
    return time_cost

#def anticipated_removal_cost(barrels, allocations):
#    costs = []
#    for barrel, allocation in zip(barrels, allocations):
#        for row, row_allocation in zip(rows, allocation):
#            costs.append(removal_time_cost(barrel, row) * row_allocation)
#    return sum(costs)
#
#def barrel_allocation_constraint(barrels, allocations):
#    return 0

#allocations = []
#for barrel in removed_barrels:
#    allocations.append([dimod.Binary(f"{barrel['label']}-{row}")
#                                    for row in range(len(rows))])
# print("allocations", allocations)

bqm = dimod.BinaryQuadraticModel(0, dimod.BINARY)

# add time costs to the model
for barrel in removed_barrels:
    for r in range(len(rows)):
        bqm.add_linear(f"{barrel['label']}-{r}", removal_time_cost(barrel, rows[r]))

print("added removal costs", bqm)

# add constraint that each removed barrel is in exactly one place
lagrange_multiplier = 3

for barrel in removed_barrels:
    for r1 in range(len(rows)):
        for r2 in range(len(rows)):
            if r1 == r2:
                bqm.add_linear(f"{barrel['label']}-{r1}", 
                               -lagrange_multiplier)
                break
            else:
                bqm.add_quadratic(f"{barrel['label']}-{r1}",
                                  f"{barrel['label']}-{r2}",
                                  2 * lagrange_multiplier)

print("added barrel uniqueness constraint", bqm)

# add constraint on row capacity using slack variables 
for r in range(len(rows)):
    spaces = row_capacity - len(rows[r])
    if not spaces or spaces >= len(removed_barrels): continue
    print(f"add constraint cap-{r} with spaces {spaces}")
    bqm.add_linear_inequality_constraint([(f"{barrel['label']}-{r}", 1) for barrel in removed_barrels],
                                         lagrange_multiplier = lagrange_multiplier,
                                         lb = 0,
                                         ub = spaces,
                                         label = f"cap-{r}")

print("added row capacity constraint", bqm)

#cqm.set_objective(anticipated_removal_cost(removed_barrels, allocations)
#                + barrel_allocation_constraint(removed_barrels, allocations))

#sampleset = EmbeddingComposite(DWaveSampler()).sample(bqm, num_reads=10)
sampleset = neal.SimulatedAnnealingSampler().sample(bqm, num_reads=10)
print(sampleset)