#%%
from pyomo.environ import *

# Create the model object
model = ConcreteModel()
size= 7

lambda_ = [0.0823295454463125,0.0823295454463125,0.0823295454463125,0.0823295454463125,0.0823295454463125,0.0823295454463125,0.0823295454463125]
U = 3
# Define the index set and decision variables
model.i = RangeSet(1, size)
model.j = RangeSet(1, size-1)

model.u = Var(model.i, initialize=0, within=NonNegativeReals)
model.x = Var(model.i,  initialize=0, within=Binary)
#%%

# Define the objective function
# model.obj = Objective(expr=sum(((model.z[j] - z_hat[1][j])**2) for i in model.i), sense=minimize)
model.obj = Objective(expr=sum((lambda_[i-1]**2)/((sum(lambda_))*2*model.u[i]*(model.u[i]-model.x[i])) for i in model.i), sense=minimize)
# z_hat[*][0] is year that is why both indices can be i


# Define the constraints

#last service rate =1
model.c1=Constraint(expr = model.u[size]==1)
#recursive service rate equation
model.c2 = ConstraintList()
for j in model.j:
    model.c2.add( model.u[j]*((model.x[j]*model.u[j+1])+((1-model.x[j])*(model.u[j+1]+1))) == model.u[j+1])
#sum of slots
model.c3 = Constraint(expr = sum(model.x[i] for i in model.i) <= U)
#steady-state condition
model.c4 = ConstraintList()
for i in model.i:
     model.c4.add(model.u[i]>=lambda_[i-1])


# Print the model
model.pprint()


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

