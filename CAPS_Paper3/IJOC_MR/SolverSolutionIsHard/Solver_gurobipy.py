#%%
import gurobipy as gp
from gurobipy import GRB

# Create a new model
m = gp.Model("mip1")
size= 7

lambda_ = [0.0823295454463125,0.0823295454463125,0.0823295454463125,0.0823295454463125,0.0823295454463125,0.0823295454463125,0.0823295454463125]
U = 3
# Define the index set and decision variables
# model.i = RangeSet(1, size)
# model.j = RangeSet(1, size-1)

# model.u = Var(model.i, initialize=0, within=NonNegativeReals)
u = m.addVars(size, lb=0.0, vtype=GRB.CONTINUOUS, name="u")
x = m.addVars(size, vtype=GRB.BINARY, name="x")


# model.x = Var(model.i,  initialize=0, within=Binary)
#%%

# Define the objective function
# model.obj = Objective(expr=sum(((model.z[j] - z_hat[1][j])**2) for i in model.i), sense=minimize)
# model.obj = Objective(expr=sum((lambda_[i-1]**2)/((sum(lambda_))*2*model.u[i]*(model.u[i]-model.x[i])) for i in model.i), sense=minimize)
# z_hat[*][0] is year that is why both indices can be i
m.setObjective(gp.quicksum ((lambda_[i-1]**2)/((sum(lambda_))*2*u[i]*(u[i]-x[i])) for i in range(size)), GRB.MINIMIZE)


# Define the constraints

#last service rate =1
m.addConstr(u[size]==1, name = 'c1')
#recursive service rate equation
m.addConstrs((u[j]*((x[j]*u[j+1])+((1-x[j])*(u[j+1]+1))) == u[j+1] for j in range(size-1)), name = 'c2')
# model.c2 = ConstraintList()
# for j in model.j:
#     model.c2.add( model.u[j]*((model.x[j]*model.u[j+1])+((1-model.x[j])*(model.u[j+1]+1))) == model.u[j+1])
#sum of slots
m.addConstr(gp.quicksum (x[i] for i in range(size)) <= U, name = 'c3')
# model.c3 = Constraint(expr = sum(model.x[i] for i in model.i) <= U)
#steady-state condition
m.addConstrs((u[i]>=lambda_[i] for i in range(size)), name = 'c4')
# model.c4 = ConstraintList()
# for i in model.i:
#      model.c4.add(model.u[i]>=lambda_[i-1] )


# Print the model
# model.pprint()


#%%
# Solve the model
solver = SolverFactory('gurobi')
solver.solve(model)
#%%
# Print the solution
best_preds=[]
for j in model.j:
    for i in model.i:
        if(model.y[i,j]()==1):
          print(f"y[{i},{j}] =  {model.y[i,j]()}")
          best_preds= best_preds + [i+1]
print(f"obj = {model.obj()}")
for j in model.j:
  print(f"z[{j}] =  {model.z[j]()}")

