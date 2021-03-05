import pyomo.environ as pe
import networkx as nx
import matplotlib.pyplot as plt
from pyomo.core.base.misc import display
from pyomo.opt.base.solvers import SolverFactory
import os

# Change to abstract 

def minlp_2bl(SD_flow=False, visualize=False):
    # PYOMO MODEL
    m = pe.ConcreteModel(name='minlp_2bl')

    # INPUTS
    S = 2  # Number of supply tanks
    B1 = 4  # Number of tanks in first blending tank
    B2 = 4  # Number of tanks in second blending tank
    D = 2  # Number of demand tanks
    T = 6  # Number of time periods
    Q = ['A']  # Substances

    # SETS
    m.S = pe.RangeSet(1, S)  # Set of supply tanks
    m.B1 = pe.RangeSet(S+1, S+B1)  # Set of first blending line tanks
    m.B2 = pe.RangeSet(S+B1+1, S+B1+B2)  # Set of second blending line tanks
    m.B = m.B1 | m.B2  # Set of blending tanks
    # Set of second blending line tanks
    m.D = pe.RangeSet(S+B1+B2+1, S+B1+B2+D)
    m.N = m.S | m.B1 | m.B2 | m.D  # Set of tanks
    m.Q = pe.Set(initialize=Q)  # Set of substances
    m.T = pe.RangeSet(1, T)  # Set of discrete time periods

    # Set of existing arcs
    if SD_flow:
        m.A = (m.S*m.B1) | (m.B1*m.B2) | (m.B2*m.D) | (m.S*m.D)
    else:
        m.A = (m.S*m.B1) | (m.B1*m.B2) | (m.B2*m.D)

    # INPUTS
    I0_n = {n: 0 for n in m.N}

    C0_qb = {(q, b): 0 for q in m.Q for b in m.B}

    IL_n = {n: 0 for n in m.N}
    IU_s = {s: 0 for s in m.S}
    IU_d = {d: 0 for d in m.D}
    IU_b1 = {b1: 30 for b1 in m.B1}
    IU_b2 = {b2: 20 for b2 in m.B2}
    IU_b = {}
    for i in (IU_b1, IU_b2):
        IU_b.update(i)
    IU_n = {}
    for i in (IU_s, IU_b, IU_d):
        IU_n.update(i)

    Fmax = 30
    FL_nn = {nn: 10 for nn in m.A}
    FU_nn = {nn: Fmax for nn in m.A}

    FDL_dt = {(11, 1): 0, (11, 2): 0, (11, 3): 15, (11, 4): 15, (11, 5): 15, (11, 6): 15,
              (12, 1): 0, (12, 2): 0, (12, 3): 15, (12, 4): 15, (12, 5): 15, (12, 6): 15}
    FDU_dt = {(d, t): Fmax for d in m.D for t in m.T}

    CL_qd = {(q, d): 0 for q in m.Q for d in m.D}
    CU_qd = {('A', 11): 0.16, ('A', 12): 1}

    CL_q = {'A': 0.06}
    CU_q = {'A': 0.26}

    CIN_qs = {('A', 1): 0.06, ('A', 2): 0.26}

    FIN_st = {(1, 1): 10, (1, 2): 10, (1, 3): 10, (1, 4): 0, (1, 5): 0, (1, 6): 0,
              (2, 1): 30, (2, 2): 30, (2, 3): 30, (2, 4): 0, (2, 5): 0, (2, 6): 0}

    betaT_s = {s: 0 for s in m.S}
    betaT_d = {11: 2, 12: 1}
    alphaN_nn = {nn: 0.1 for nn in m.A}
    betaN_nn = {nn: 0 for nn in m.A}

    big_MC = 1

    # PARAMETERS

    # Initial inventories
    m.I0_n = pe.Param(m.N, initialize=I0_n, within=pe.NonNegativeReals)

    # Initial composition
    m.C0_qb = pe.Param(m.Q, m.B, initialize=C0_qb, within=pe.NonNegativeReals)

    # Inventory bounds
    m.IL_n = pe.Param(m.N, initialize=IL_n, within=pe.NonNegativeReals)
    m.IU_n = pe.Param(m.N, initialize=IU_n, within=pe.NonNegativeReals)

    # Flow bounds
    m.FL_nn = pe.Param(m.A, initialize=FL_nn, within=pe.NonNegativeReals)
    m.FU_nn = pe.Param(m.A, initialize=FU_nn, within=pe.NonNegativeReals)

    # Demanded flow bounds
    m.FDL_dt = pe.Param(m.D, m.T, initialize=FDL_dt,
                        within=pe.NonNegativeReals)
    m.FDU_dt = pe.Param(m.D, m.T, initialize=FDU_dt,
                        within=pe.NonNegativeReals)

    # Demanded composition bounds
    m.CL_qd = pe.Param(m.Q, m.D, initialize=CL_qd, within=pe.NonNegativeReals)
    m.CU_qd = pe.Param(m.Q, m.D, initialize=CU_qd, within=pe.NonNegativeReals)

    # Composition bounds
    m.CL_q = pe.Param(m.Q, initialize=CL_q, within=pe.NonNegativeReals)
    m.CU_q = pe.Param(m.Q, initialize=CU_q, within=pe.NonNegativeReals)

    # Supply conditions
    m.CIN_qs = pe.Param(m.Q, m.S, initialize=CIN_qs,
                        within=pe.NonNegativeReals)
    m.FIN_st = pe.Param(m.S, m.T, initialize=FIN_st,
                        within=pe.NonNegativeReals)

    # Economic parameters
    m.betaT_s = pe.Param(m.S, initialize=betaT_s, within=pe.NonNegativeReals)
    m.betaT_d = pe.Param(m.D, initialize=betaT_d, within=pe.NonNegativeReals)
    m.alphaN_nn = pe.Param(m.A, initialize=alphaN_nn,
                           within=pe.NonNegativeReals)
    m.betaN_nn = pe.Param(m.A, initialize=betaN_nn, within=pe.NonNegativeReals)

    # CONTINUOUS VARIABLES
    m.F = pe.Var(m.A, m.T, within=pe.NonNegativeReals)
    m.FD = pe.Var(m.D, m.T, within=pe.NonNegativeReals)
    m.I = pe.Var(m.N, m.T, within=pe.NonNegativeReals)
    m.C = pe.Var(m.Q, m.B, m.T, within=pe.NonNegativeReals)

    # BINARY VARIABLES
    m.X = pe.Var(m.A, m.T, within=pe.Binary)

    # CONSTRAINTS

    # Flow activation
    m.flow_activation = pe.ConstraintList()
    for nn in m.A:
        for t in m.T:
            m.flow_activation.add(m.F[nn, t] <= m.FU_nn[nn] * m.X[nn, t])
            m.flow_activation.add(m.F[nn, t] >= m.FL_nn[nn] * m.X[nn, t])

    # Satisfy specifications
    @m.Constraint(m.Q, m.D, m.B2, m.T)
    def satisfy_specs_1(m, q, d, b, t):
        if t > 1:
            return m.C[q, b, t-1] <= m.CU_qd[q, d] + big_MC * (1 - m.X[(b, d), t])
        else:
            return pe.Constraint.Skip

    @m.Constraint(m.Q, m.D, m.B2, m.T)
    def satisfy_specs_2(m, q, d, b, t):
        if t > 1:
            return m.C[q, b, t-1] >= m.CL_qd[q, d] - big_MC * (1 - m.X[(b, d), t])
        else:
            return pe.Constraint.Skip

    # Satisfy SD specifications
    if SD_flow:
        @m.Constraint(m.Q, m.D, m.S, m.T)
        def sd_specs_1(m, q, d, s, t):
            return m.CIN_qs[q, s] <= m.CU_qd[q, d] + big_MC * (1 - m.X[(s, d), t])

        @m.Constraint(m.Q, m.D, m.S, m.T)
        def sd_specs_2(m, q, d, s, t):
            return m.CIN_qs[q, s] >= m.CL_qd[q, d] - big_MC * (1 - m.X[(s, d), t])

    # Supply inventory balance
    if SD_flow:
        @m.Constraint(m.S, m.T)
        def supply_bal_sd(m, s, t):
            if t == 1:
                return m.I[s, t] == m.I0_n[s] + m.FIN_st[s, t] - sum(m.F[(s, n), t] for n in m.B1) - sum(m.F[(s, d), t] for d in m.D)
            else:
                return m.I[s, t] == m.I[s, t-1] + m.FIN_st[s, t] - sum(m.F[(s, n), t] for n in m.B1) - sum(m.F[(s, d), t] for d in m.D)
    else:
        @m.Constraint(m.S, m.T)
        def supply_bal(m, s, t):
            if t == 1:
                return m.I[s, t] == m.I0_n[s] + m.FIN_st[s, t] - sum(m.F[(s, n), t] for n in m.B1)
            else:
                return m.I[s, t] == m.I[s, t-1] + m.FIN_st[s, t] - sum(m.F[(s, n), t] for n in m.B1)

    # Blending inventory balances
    @m.Constraint(m.B1, m.T)
    def blend_1_bal(m, b, t):
        if t == 1:
            return m.I[b, t] == m.I0_n[b] + sum(m.F[(s, b), t] for s in m.S) - sum(m.F[(b, n), t] for n in m.B2)
        else:
            return m.I[b, t] == m.I[b, t-1] + sum(m.F[(s, b), t] for s in m.S) - sum(m.F[(b, n), t] for n in m.B2)

    @m.Constraint(m.B2, m.T)
    def blend_2_bal(m, b, t):
        if t == 1:
            return m.I[b, t] == m.I0_n[b] + sum(m.F[(n, b), t] for n in m.B1) - sum(m.F[(b, d), t] for d in m.D)
        else:
            return m.I[b, t] == m.I[b, t-1] + sum(m.F[(n, b), t] for n in m.B1) - sum(m.F[(b, d), t] for d in m.D)

    # Demand inventory balance
    if SD_flow:
        @m.Constraint(m.D, m.T)
        def demand_bal_sd(m, d, t):
            if t == 1:
                return m.I[d, t] == m.I0_n[d] + sum(m.F[(n, d), t] for n in m.B2) + sum(m.F[(s, d), t] for s in m.S) - m.FD[d, t]
            else:
                return m.I[d, t] == m.I[d, t-1] + sum(m.F[(n, d), t] for n in m.B2) + sum(m.F[(s, d), t] for s in m.S) - m.FD[d, t]
    else:
        @m.Constraint(m.D, m.T)
        def demand_bal(m, d, t):
            if t == 1:
                return m.I[d, t] == m.I0_n[d] + sum(m.F[(n, d), t] for n in m.B2) - m.FD[d, t]
            else:
                return m.I[d, t] == m.I[d, t-1] + sum(m.F[(n, d), t] for n in m.B2) - m.FD[d, t]

    # Variable implications
    @m.Constraint(m.B1, m.S, m.B2, m.T)
    def implication_b1(m, b1, s, b2, t):
        return m.X[(s, b1), t] + m.X[(b1, b2), t] <= 1

    @m.Constraint(m.B1, m.D, m.B2, m.T)
    def implication_b2(m, b1, d, b2, t):
        return m.X[(b1, b2), t] + m.X[(b2, d), t] <= 1

    # Bounds
    m.I_bounds = pe.ConstraintList()
    for n in m.N:
        for t in m.T:
            m.I_bounds.add(m.IL_n[n] <= m.I[n, t])
            m.I_bounds.add(m.I[n, t] <= m.IU_n[n])

    m.FD_bounds = pe.ConstraintList()
    for d in m.D:
        for t in m.T:
            m.I_bounds.add(m.FDL_dt[d, t] <= m.FD[d, t])
            m.I_bounds.add(m.FD[d, t] <= m.FDU_dt[d, t])

    m.C_bounds = pe.ConstraintList()
    for b in m.B:
        for q in m.Q:
            for t in m.T:
                m.C_bounds.add(m.CL_q[q] <= m.C[q, b, t])
                m.C_bounds.add(m.C[q, b, t] <= m.CU_q[q])

    # OBJECTIVE
    if SD_flow:
        @m.Objective(sense=pe.maximize)
        def obj(m):
            return sum(sum(m.betaT_d[d] * m.F[(n, d), t] for d in m.D for n in m.B2) + sum(m.betaT_d[d] * m.F[(s, d), t] for d in m.D for s in m.S) - sum(m.alphaN_nn[nn] * m.X[nn, t] + m.betaN_nn[nn] * m.F[nn, t] for nn in m.A) for t in m.T)
    else:
        @m.Objective(sense=pe.maximize)
        def obj(m):
            return sum(sum(m.betaT_d[d] * m.F[(n, d), t] for d in m.D for n in m.B2) - sum(m.alphaN_nn[nn] * m.X[nn, t] + m.betaN_nn[nn] * m.F[nn, t] for nn in m.A) for t in m.T)

    return m


def solver(m):
    dir_path = os.path.dirname(os.path.abspath(__file__))
    gams_path = os.path.join(dir_path, "gamsfiles/")
    if not(os.path.exists(gams_path)):
        print('Directory for automatically generated files ' +
              gams_path + ' does not exist. We will create it')
        os.makedirs(gams_path)

    # SOLVE
    solvername = 'gams'
    opt = SolverFactory(solvername, solver='cplex')
    results = opt.solve(m, tee=True,
                        # Uncomment the following lines if you want to save GAMS models
                        # keepfiles=True,
                        # tmpdir=gams_path,
                        # symbolic_solver_labels=True,

                        add_options=[
                            'option reslim = 200;'
                            'option optcr = 0.01;'
                            # Uncomment the following lines to setup IIS computation of BARON through option file
                            # 'GAMS_MODEL.optfile = 1;'
                            # '\n'
                            # '$onecho > baron.opt \n'
                            # 'CompIIS 1 \n'
                            # '$offecho'
                        ])

    print('Objective:', round(pe.value(m.obj), 5))
    return m


def visualize(m, SD_flow=False):
    x = [1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4]
    y = [3, 2, 4, 3, 2, 1, 4, 3, 2, 1, 3, 2]
    pos = {(i+1): (x[i], y[i]) for i in range(len(m.N))}

    for t in m.T:
        print('Optimal operating flows at the end of period', t, 'are:')
        arcs = []
        flows = []
        for nn in m.A:
            if pe.value(m.F[nn, t]) >= 0.01:
                arcs.append(((nn)))
                flows.append(round(pe.value(m.F[nn, t]), 2))

        graph = nx.DiGraph()
        for i in range(len(arcs)):
            graph.add_edge(arcs[i][0], arcs[i][1])

        pairs = list(zip(list(arcs), list(flows)))
        edgelabels = dict(pairs)

        nodelabels = {i: str(i) for i in range(1, len(m.N)+1)}

        nx.draw_networkx(graph, pos, node_size=700, node_color='skyblue',
                         width=1.5, nodelist=list(range(1, len(m.N)+1)), with_labels=True)
        nx.draw_networkx_edge_labels(graph, pos, edge_labels=edgelabels)
        nx.draw_networkx_labels(graph, pos, nodelabels)
        plt.show()
        if SD_flow:
            for s in m.S:
                for d in m.D:
                    if pe.value(m.F[(s, d), t]) >= 0.1:
                        for q in m.Q:
                            print('Specification of', q, 'that exited through', d, 'from', s, 'is', round(
                                m.CIN_qs[q, s], 3), 'and it must be between [', m.CL_qd[q, d], ',', m.CU_qd[q, d], ']')

        if t > 1:
            for d in m.D:
                for b2 in m.B2:
                    if pe.value(m.F[(b2, d), t]) >= 0.1:
                        for q in m.Q:
                            print('Specification of', q, 'that exited through', d, 'from', b2, 'is', round(pe.value(
                                m.C[q, b2, t-1], 3)), 'and it must be between [', m.CL_qd[q, d], ',', m.CU_qd[q, d], ']')
        print()


if __name__ == "__main__":
    m = minlp_2bl(SD_flow=False)
    m_solved = solver(m)
    visualize(m_solved, SD_flow=False)
