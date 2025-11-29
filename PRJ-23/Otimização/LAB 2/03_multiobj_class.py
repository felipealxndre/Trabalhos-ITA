'''
INSTITUTO TECNOLÓGICO DE AERONÁUTICA
PROGRAMA DE ESPECIALIZAÇÃO EM ENGENHARIA AERONÁUTICA
OTIMIZAÇÃO MULTIDISCIPLINAR

Código com aplicações de ferramentas do Python para otimização
multi-objetivo (pymoo 0.6.0)

Maj. Ney Sêcco 2025
'''

# IMPORTS
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.core.problem import ElementwiseProblem
import numpy as np
import matplotlib.pyplot as plt

# EXECUTION

# Define optimization problem
class MyProblem(ElementwiseProblem):

    def __init__(self):

        # Set general characteristics of the problem
        super().__init__(n_var=2, # Number of design variables - número de variáveis de projeto, x1 e x2
                         n_obj=2, # Number of objective functions - número de objetivos, f1 e f2
                         n_constr=2, # Number of constraints
                         # batentes
                         xl=[-2, -2], # Lower bounds of design variables
                         xu=[2, 2]) # Upper bound of design variables

    def _evaluate(self, X, out, *args, **kwargs):

        # Split design variables
        x1 = X[0]  # x1 é o primeiro elemento de X
        x2 = X[1]  # x2 é o segundo elemento de X

        # Compute objectives
        f1 = x1**2 + x2**2
        f2 = (x1-1)**2 + x2**2

        # Compute constraints - aqui tem que ser menor igual a zero
        g1 = 2*(x1-0.1) * (x1-0.9) / 0.18
        g2 = - 20*(x1-0.4) * (x1-0.6) / 4.8

        # Gather results
        out["F"] = [f1, f2]
        out["G"] = [g1, g2]

# Create an instance of the Problem
problem = MyProblem()

# Select optimization algorithm
# tamanho da população: pop size
algorithm = NSGA2(pop_size=100, eliminate_duplicates=True)

# Solve the optimization
res = minimize(problem,
               algorithm,
               ('n_gen', 100), # número de gerações - quantas vezes ele vai replicar
               seed=1,
               verbose=True)


# Plot the Pareto in both objective and design spaces
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Plot the Pareto in the objective space
ax1.plot(res.F[:,0], res.F[:,1], 'o')
ax1.set_title('Objective Space', fontsize=20)
ax1.set_xlabel(r'$f_1$', fontsize=20)
ax1.set_ylabel(r'$f_2$', fontsize=20)

# Plot the Pareto in the design space
ax2.plot(res.X[:,0], res.X[:,1], 'o')
ax2.set_title('Design Space', fontsize=20)
ax2.set_xlabel(r'$x_1$', fontsize=20)
ax2.set_ylabel(r'$x_2$', fontsize=20)
ax2.axis('equal')

plt.tight_layout()
plt.show()
#fig.savefig('pareto_nsga.pdf')
