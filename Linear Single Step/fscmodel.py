# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 13:49:17 2018
@author: j.robinson
"""

from __future__ import division
from pyomo.environ import *
from pyomo.opt import SolverFactory
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

input_file = 'input_BAT2030.xlsx'

class Source:
    def __init__(self, name, energyType, capex, opex, CO2, minProd, maxProd, pos):
        self.name = name
        self.energyType = energyType
        self.capex = capex
        self.opex = opex
        self.CO2 = CO2
        self.outcons = []
        self.minProd = minProd
        self.maxProd = maxProd
        self.pos = pos
    
    def __str__(self):
        return "Source:" + self.name + ", " + self.energyType
    
    def __lt__(self,other):
        if isinstance(other, Source):
            return self.name < other.name


class Sink:
    def __init__(self, name, capex, opex, energyType, demand, pos):
        self.name = name
        self.energyType = energyType
        self.capex = capex
        self.opex = opex
        self.demand = demand
        self.incons = []
        self.pos = pos
        
    def __str__(self):
        return "Sink:" + self.name + ", " + self.energyType
    
    def __lt__(self, other):
        if isinstance(other, Sink):
            return self.name < other.name
    
class Transformer:
    def __init__(self, name, capex, opex, totalEff, outMin, outMax, pos):
        self.name = name
        self.capex = capex
        self.opex = opex
        self.totalEff = totalEff
        self.outMin = outMin
        self.outMax = outMax
        self.inputs = {}
        self.products = {}
        self.incons = []
        self.outcons = []
        self.pos = pos
    
    def __str__(self):
        return "Transformer:" + self.name
    
    def __lt__(self, other):
        if isinstance(other, Transformer):
            return self.name < other.name
        
class Hub:
    def __init__(self, name, energyType, capex, opex, pos):
        self.name = name
        self.energyType = energyType
        self.capex = capex
        self.opex = opex
        self.incons = []
        self.outcons = []
        self.pos = pos
    
    def __str__(self):
        return "Hub:" + self.name + ", " + self.energyType
    
    def __lt__(self, other):
        if isinstance(other, Hub):
            return self.name < other.name

class Connection:
    def __init__(self, name, inp, out, energyType):
        self.name = name
        self.inp = inp
        self.out = out
        self.energyType = energyType
    
    def __lt__(self, other):
        if isinstance(other, Connection):
            return self.name < other.name
    def __str__(self):
        return "Connection:" + self.name + ", " + self.energyType


def createModel(SourceList, SinkList, TransList, ConnList, HubList, CO2):
    M = ConcreteModel()
    
    M.connectors = Set(initialize = ConnList)
    M.sources = Set(initialize = SourceList)
    M.sinks = Set(initialize = SinkList)
    M.trans = Set(initialize = TransList)
    M.hubs = Set(initialize = HubList)
    M.stations = Set(initialize = SourceList + SinkList + TransList + HubList)
    
    M.c = Param(M.stations, mutable = True)
    M.carbon = Param(M.sources, mutable = True)
    M.cape = Param(M.stations, mutable = True)
    
    # For the amount in facilities, for calculating Opex. For transformer, the amount coming out
    M.facilities = Var( M.stations, domain = NonNegativeReals)
    # Whether a facility is being used. For calculating Capex
    M.isopen = Var(M.stations, domain = Boolean)
    # Amount going through connectors
    M.connections = Var(M.connectors, domain = NonNegativeReals)
    # Amount coming into a transformer. Used to consider transformer production ratio
    M.trintotals = Var(M.trans, domain = NonNegativeReals)
    
    # Populates capex costs
    for fac in M.stations:
        M.cape[fac]=fac.capex
    
    # Constructs cost vector from opex and carbon constraints from sources.
    for fac in M.stations:
        M.c[fac] = fac.opex
        if isinstance(fac, Source):
            M.carbon[fac] = fac.CO2
    
    
    def sourcecount(model, source):
        return M.facilities[source] == sum(M.connections[con] for con in source.outcons)
    
    M.sourcesum = Constraint(M.sources, rule = sourcecount)
    
    for source in M.sources:
        M.facilities[source].setub(source.maxProd)
        M.facilities[source].setlb(source.minProd)
    
    def transrule(model, tra):
        return M.facilities[tra] == tra.totalEff * M.trintotals[tra]
    
    def transcount(model, tra):
        return M.trintotals[tra] ==  sum(M.connections[con] for con in tra.incons)
    
    def inputratiorule(model, con):
        for tra in TransList:
            if con in tra.incons:
                return tra.inputs[con.energyType] * M.trintotals[tra] == M.connections[con]
        return Constraint.Skip
     
    def productratiorule(model, con):
        for tra in TransList:
            if con in tra.outcons:
                etype = con.energyType
                return tra.products[etype] * M.facilities[tra] == M.connections[con]
        return Constraint.Skip
    
    for tran in M.trans:
        M.facilities[tran].setub(tran.outMax)
        M.facilities[tran].setlb(tran.outMin)

    M.transconstraint = Constraint(M.trans, rule = transrule)
    M.transsum = Constraint(M.trans, rule = transcount)
    M.inputconstraint = Constraint(M.connectors, rule = inputratiorule)
    M.productconstraint = Constraint(M.connectors, rule = productratiorule)
    
    
    def sinkrule(model, sink):
        return sum(M.connections[con] for con in sink.incons) == sink.demand
    
    def sinkcount(model,sink):
        return M.facilities[sink]== sum(M.connections[con] for con in sink.incons)
    
    M.sinkconstraint = Constraint(M.sinks, rule = sinkrule)
    M.sinksum = Constraint(M.sinks, rule = sinkcount)
    
    
    def hubrule(model, hub):
        return sum(M.connections[con] for con in hub.incons)==sum(M.connections[con] for con in hub.outcons)
    
    def hubcount(model,hub):
        return M.facilities[hub] == sum(M.connections[con] for con in hub.incons)
    
    M.hubconstraint = Constraint(M.hubs, rule = hubrule)
    M.hubsum = Constraint(M.hubs, rule = hubcount)
    
    # Quadratic constraint that turns isopen on and off
    def binrule(model, fac):
        return M.facilities[fac] - M.isopen[fac]*M.facilities[fac] <= 0
    
    M.checkopen = Constraint(M.stations, rule = binrule)


    M.Co2limit = Constraint(expr = summation(M.facilities,M.carbon,index = M.sources) <= CO2)    
        
    def objrule(model):
       ob = summation(model.facilities,model.c, index=M.stations) + summation(model.cape, model.isopen, index=M.stations)
       return ob
            
    M.Obj = Objective(rule = objrule, sense = minimize)
            
    return M

def opti(model):
    opt = SolverFactory('gurobi')
    results = opt.solve(model, tee = True)
    print(model.display())
    return results


def checkModel(ConnList, entypes):
    for con in ConnList:
        if con.energyType not in entypes:
            raise ValueError(str(con) + ' has an unrecognized energy type.')
    
    for Source in SourceList:
       if not Source.outcons:
           print('\nWARNING: ' + Source.name + ' has empty out connections, so it probably is not being used. Would you like to check that?' + '\n')   
    #What more can be added?
    return None

SourceIn    = pd.read_excel(input_file, 'Sources', index_col=None, na_values=['NA'])
SinkIn      = pd.read_excel(input_file, 'Sinks', index_col=None, na_values=['NA'])
TransIn     = pd.read_excel(input_file, 'Transformers', index_col=None, na_values=['NA'])
HubIn       = pd.read_excel(input_file, 'Hubs', index_col=None, na_values=['NA'])
ConnIn      = pd.read_excel(input_file, 'Connectors', index_col=None, na_values=['NA'])
RestrIn     = pd.read_excel(input_file, 'Restrictions', index_col=None, na_values=['NA'])

SourceList     = []
SinkList       = []
TransList      = []
HubList        = []
ConnList       = []
FuelTypeList   = []
DemandTypeList = []
outcols = ['Total Cost']


# Import restrictions, just CO2 for now
CO2Max = RestrIn.loc[0,'CO2 Max']

# Energy sources available from sources
for i in range(len(SourceIn.index)):
    if not SourceIn.loc[i,'EnergyType'] in FuelTypeList:
        FuelTypeList.append(SourceIn.loc[i,'EnergyType'])
        outcols.append(SourceIn.loc[i,'EnergyType'])
        
# Energy types demanded at sinks     
for i in range(len(SinkIn.index)):
    if not SinkIn.loc[i, 'EnergyType'] in DemandTypeList:
        DemandTypeList.append(SinkIn.loc[i, 'EnergyType'])

# All energy types 
EnergyList = FuelTypeList + DemandTypeList

G = nx.DiGraph()
posits = {}
labelpos = {}

# Initialize the connectors        
for i in range(len(ConnIn.index)):
    ConnList.append(Connection(name = ConnIn.loc[i,'Name'],
                              inp = ConnIn.loc[i,'From'],
                              out = ConnIn.loc[i,'To'],
                              energyType = ConnIn.loc[i,'EnergyType']))

# Initialize the Sources
for i in range(len(SourceIn.index)):
    SourceList.append(Source(name = SourceIn.loc[i,'Name'],
                             energyType = SourceIn.loc[i,'EnergyType'],
                             capex=SourceIn.loc[i,'Capex'], 
                             opex = SourceIn.loc[i,'Opex'], 
                             CO2 = SourceIn.loc[i,'CO2'],
                             minProd = SourceIn.loc[i,'MinProduction'],
                             maxProd = SourceIn.loc[i, 'MaxProduction'],
                             pos = (SourceIn.loc[i, 'X'],SourceIn.loc[i, 'Y'])))
    G.add_node(SourceList[i].name, pos = SourceList[i].pos, shape = 's', color = 'g')
    posits[SourceList[i].name] =  SourceList[i].pos
    labelpos[SourceList[i].name] = ((SourceIn.loc[i, 'X'],SourceIn.loc[i, 'Y'] - 10))
    
    for con in ConnList:
        if con.inp==SourceList[i].name and con.energyType==SourceList[i].energyType:
            SourceList[i].outcons.append(con)

# Initialize the sinks
for i in range(len(SinkIn.index)):
    SinkList.append(Sink(name = SinkIn.loc[i,'Name'],
                         energyType = SinkIn.loc[i,'EnergyType'],
                         capex = SinkIn.loc[i,'Capex'],
                         opex = SinkIn.loc[i,'Opex'],
                         demand = SinkIn.loc[i,'Demand'],
                         pos = (SinkIn.loc[i, 'X'],SinkIn.loc[i, 'Y'])))
    G.add_node(SinkList[i].name, pos= SinkList[i].pos, shape = 's', color = 'r')
    posits[SinkList[i].name] =  SinkList[i].pos
    labelpos[SinkList[i].name] = (SinkIn.loc[i, 'X'],SinkIn.loc[i, 'Y'] - 10)
    
    for con in ConnList:
        if con.out==SinkList[i].name and con.energyType==SinkList[i].energyType:
            SinkList[i].incons.append(con)

# Initialize the transfomers
for i in range(len(TransIn.index)):
    TransList.append(Transformer(name = TransIn.loc[i,'Name'],
                                 capex = TransIn.loc[i,'Capex'],
                                 opex = TransIn.loc[i,'Opex'],
                                 totalEff = TransIn.loc[i,'TotalEff'],
                                 outMin = TransIn.loc[i, 'OutMin'],
                                 outMax = TransIn.loc[i, 'OutMax'],
                                 pos = (TransIn.loc[i, 'X'],TransIn.loc[i, 'Y'])))
    G.add_node(TransList[i].name, pos = TransList[i].pos, shape = '8', color = 'b')
    posits[TransList[i].name] =  TransList[i].pos
    labelpos[TransList[i].name] = (TransIn.loc[i, 'X'],TransIn.loc[i, 'Y'] - 10)
    
    
    outcols.append(TransList[i].name + 'Production')
    
    k = 0
    
    for j in range(len(TransIn.loc[i,'Input0':'Prod0'])-1):
        x = int(j/2)
        inp = TransIn.loc[i,'Input'+str(x)]       
        if k % 2 == 0 and isinstance(inp,str):
            if not inp in EnergyList:
                EnergyList.append(inp)
            TransList[i].inputs[inp] = TransIn.loc[i,'InRatio'+str(x)]
        k = k + 1
        
    k = 0
    
    for j in range(len(TransIn.loc[i,'Prod0':])):
        x = int(j/2)
        product = TransIn.loc[i,'Prod'+str(x)]       
        if k % 2 == 0 and isinstance(product,str):
            if not product in EnergyList:
                EnergyList.append(product)
            TransList[i].products[product] = TransIn.loc[i,'SubEff'+str(x)]
            x = x + 1
        k = k + 1
 
    for con in ConnList:
        if con.out==TransList[i].name and con.energyType in TransList[i].inputs:
            TransList[i].incons.append(con)
        elif con.inp==TransList[i].name and con.energyType in TransList[i].products:
            TransList[i].outcons.append(con)
            outcols.append(TransList[i].name + '-' + con.energyType)

# Initialize the Hubs   
for i in range(len(HubIn.index)):
    HubList.append(Hub(name = HubIn.loc[i,'Name'],
                       energyType = HubIn.loc[i,'EnergyType'],
                       capex = HubIn.loc[i,'Capex'],
                       opex = HubIn.loc[i,'Opex'],
                       pos = (HubIn.loc[i, 'X'],HubIn.loc[i, 'Y'])))
    G.add_node(HubList[i].name, pos = HubList[i].pos, shape = 'o', color = 'white')
    posits[HubList[i].name] =  HubList[i].pos
    labelpos[HubList[i].name] = (HubIn.loc[i, 'X'],HubIn.loc[i, 'Y'] - 10)
    
    outcols.append(HubList[i].name + 'Usage')
    for con in ConnList:
        if con.out==HubList[i].name and con.energyType==HubList[i].energyType:
            HubList[i].incons.append(con)
        elif con.inp==HubList[i].name and con.energyType==HubList[i].energyType:
            HubList[i].outcons.append(con)
    


checkModel(ConnList, EnergyList)

model = createModel(SourceList, SinkList, TransList, ConnList, HubList, CO2 = CO2Max)

results = opti(model)


# Output formatting starts here

    
outdf = pd.DataFrame(np.zeros((1,len(outcols))), columns = outcols)
outdf.at[0, 'Total Cost'] = model.Obj()
for fac in model.stations:
    if isinstance(fac, Source):
        outdf.at[0, fac.energyType] = model.facilities[fac].value
    elif isinstance(fac, Transformer):
        outdf.at[0, fac.name + 'Production'] = model.facilities[fac].value
        for con in fac.outcons:
            outdf.at[0, fac.name + '-' + con.energyType] = model.connections[con].value
    else:
        outdf.at[0, fac.name + 'Usage'] = model.facilities[fac].value
        


outdf.to_excel('output.xlsx', sheet_name='Sheet1')

plt.figure(figsize = (15,9))
plt.axis('off') 
for con in ConnList:
    if model.connections[con].value > 0.000001:
        G.add_edge(con.inp, con.out)
        nx.draw_networkx_edges(G, pos = posits, edgelist = [(con.inp, con.out)], width = model.connections[con].value/60)
    else:
        G.add_edge(con.inp, con.out)

nx.draw_networkx_labels(G, pos = labelpos)
for node in G.nodes():
    nx.draw_networkx_nodes(G, pos = posits, nodelist = [node], node_color = G.node[node]['color'], node_shape = G.node[node]['shape'])

plt.show()
plt.figure(figsize = (9, 5))
nx.draw(G, pos = posits, with_labels = True, node_shape = 's')
checkModel(ConnList, EnergyList)
