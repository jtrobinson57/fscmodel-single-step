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
import math

class Source:
    def __init__(self,name,energyType,capex,opex,CO2,minProd,maxProd):
        self.name = name
        self.energyType = energyType
        self.capex = capex
        self.opex = opex
        self.CO2 = CO2
        self.outcons = []
        self.minProd = minProd
        self.maxProd = maxProd
    
    def __str__(self):
        return "Source:" + self.name + ", " + self.energyType
    
    def __lt__(self,other):
        if isinstance(other, Source):
            return self.name < other.name


class Sink:
    def __init__(self,name,capex,opex,energyType,demand):
        self.name = name
        self.energyType = energyType
        self.capex = capex
        self.opex = opex
        self.demand = demand
        self.incons = []
        
    def __str__(self):
        return "Sink:" + self.name + ", " + self.energyType
    
    def __lt__(self,other):
        if isinstance(other, Sink):
            return self.name < other.name
    
class Transformer:
    def __init__(self, name, capex, opex, totalEff, outMin, outMax, CO2Ratio):
        self.name = name
        self.capex = capex
        self.opex = opex
        self.totalEff = totalEff
        self.outMin = outMin
        self.outMax = outMax
        self.CO2Ratio = CO2Ratio
        self.inputs = {}
        self.products = {}
        self.incons = []
        self.outcons = []
        self.hynum = None
    
    def __str__(self):
        return "Transformer:" + self.name
    
    def __lt__(self,other):
        if isinstance(other, Transformer):
            return self.name < other.name
        
class Hub:
    def __init__(self,name,energyType,capex=0,opex=0):
        self.name = name
        self.energyType = energyType
        self.capex = capex
        self.opex = opex
        self.incons = []
        self.outcons = []
    
    def __str__(self):
        return "Hub:" + self.name + ", " + self.energyType
    
    def __lt__(self,other):
        if isinstance(other, Hub):
            return self.name < other.name

class Connection:
    def __init__(self,name,inp,out,energyType):
        self.name = name
        self.inp = inp
        self.out = out
        self.energyType = energyType
    
    def __lt__(self,other):
        if isinstance(other, Connection):
            return self.name < other.name
    def __str__(self):
        return "Connection:" + self.name + ", " + self.energyType

class CO2Loc:
    def __init__(self,name,ind,postal,dist,cap,costPKG):
        self.name = name
        self.ind = ind
        self.postal = postal
        self.dist = dist / 100.0       #Hundreds of km
        self.cap = cap #/ 1000           #MW
        self.capPJ = 0
        self.capex = 0                 #Euros
        self.indOpex = 0               #Euros
        self.dirOpex = 0               #Euros per kg H2
        self.K = 0                     #Euros
        self.costPKG = costPKG
        
    def findCapex(self):
        
        self.capex = 665 - (349.721)*(1-math.e**(-0.015056*self.cap / 1000))    #returns euro/KW
        self.capex = self.capex * 1000 * self.cap                        #returns euros
        
    def findDirOpex(self):
        
        self.dirOpex = 5.190866 + (3.999796 - 5.190866)/(1 + (self.dist/2.020612)**1.534203) #returns euros/kg of H2
        self.dirOpex = self.dirOpex / 120 * 10**9          #converts to euros/PJ of fuel
        
    def findIndOpex(self):
        
        #I only made this so verbose to make the unit conversions a bit more clear
        
        MW = self.cap# * 1000
        MJpa = MW * 3600 * 8000  #Converted to MJ/a
        MJph = MJpa / 8000       #Converted to MJ/h
        KGph = MJph / 43.1       #Converted to KG/h
        
        self.indOpex = 2.13 * ((KGph)**0.242) * 4 * (8000/24) * 37.32 * 4  #returns euros 
    
    def checkMinMax(self):     #Checks to make sure plant capacity is within given bounds
                               #according to the input restrictions sheet       
        if self.cap > maxCO2PlantSize:     #Admittedly, this function only checks maxes, mins are checked below
            self.cap = maxCO2PlantSize     #at the point of read-in
    
    def changeCapUnit(self):
        
        self.capPJ = self.cap * 3600 * 8000 / 1000000000  #Switch from MW to PJ/yr for a single year
    
    def __lt__(self,other):
        if isinstance(other, Connection):
            return self.name < other.name
    def __str__(self):
        return "CO2 Location:" + self.name
    

    
def createModel(SourceList, SinkList, TransList, ConnList, HubList, CO2LocList, CO2):

    M = ConcreteModel()
    
    M.connectors = Set(initialize = ConnList)
    M.sources = Set(initialize = SourceList)
    M.sinks = Set(initialize = SinkList)
    M.trans = Set(initialize = TransList)
    M.hytrans = Set(initialize = H2TransList)
    M.hubs = Set(initialize = HubList)
    M.stations = Set(initialize = SourceList + SinkList + TransList + HubList)
    M.locations = Set(initialize = CO2LocList)

    
    
    letsinit = {}
    for loc in M.locations:
        iii = []
        for j in range(hyn):
            iii.append(locationNum*j + loc.ind)
        letsinit[loc] = iii
    M.letsgo = Set(M.locations, initialize = letsinit)
#    print(M.letsgo[CO2LocList[1]].value)
    
    M.c = Param(M.stations, mutable = True)
    M.carbon = Param(M.sources, mutable = True)
    M.cape = Param(M.stations, mutable = True)
    M.loccap = Param(M.locations, mutable = True)
    M.locopex = Param(M.locations, mutable = True)
    M.CO2costs = Param(M.locations, mutable = True)
    
    #For the amount in facilities, for calculating Opex. For transformer, the amount coming out
    M.facilities = Var( M.stations, domain = NonNegativeReals)
    #Whether a facility is being used. For calculating Capex
    M.isopen = Var(M.stations, domain = Boolean)
    #Amount going through connectors
    M.connections = Var(M.connectors, domain = NonNegativeReals)
    #Amount coming into a transformer. Used to consider transformer production ratio
    M.trintotals = Var(M.trans, domain = NonNegativeReals)
    M.carbonsum = Var(domain = NonNegativeReals)
    M.carbonsum.setub(CO2)
    #Amount of hydrogen in MJ going into C02 locs.
    M.hydrouse = Var(M.locations, domain = NonNegativeReals)
    #Whether the CO2 locations are open
    M.locopen = Var(M.locations, domain = Boolean)
    #Assignment of C02locs to H2 transformers
    M.assignments = Var(range(locationNum*len(H2TransList)), domain = Boolean)
    
    #Populates capex costs
    for fac in M.stations:
        M.cape[fac]=fac.capex
        
    for loc in M.locations:
        M.loccap[loc] = loc.K + (loc.indOpex)    #No longer using the MW cap
        M.locopex[loc] = loc.dirOpex 
        
    #Constructs cost vector from opex and carbon constraints from sources.
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
        
    #Transformer related equations.
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
#Dealing with CO2 locations here    
#Set maximum.
    for loc in CO2LocList:
        M.hydrouse[loc].setub(loc.capPJ)
    
#Sum equation that adds up numbers **to the hydrogen In connector.
    def hydrosum(model, hy):
        for con in hy.incons:
            if con.energyType=='hydrogen':
                return sum(model.assignments[hy.hynum*locationNum + loc.ind] * model.hydrouse[loc] for loc in CO2LocList) >= M.connections[con]
    
    M.hopethisworks = Constraint(M.hytrans, rule = hydrosum)

    def binruleloc(model, loc):
        return model.hydrouse[loc] - model.locopen[loc]*model.hydrouse[loc] <= 0
    
    M.checklocopen = Constraint(M.locations, rule = binruleloc)

    M.sososo = SOSConstraint(M.locations, var = M.assignments, index = M.letsgo, sos = 1 )
    
    #Quadratic constraint that turns isopen on and off
    def binrule(model, fac):
        return M.facilities[fac] - M.isopen[fac]*M.facilities[fac] <= 0
    
    M.checkopen = Constraint(M.stations, rule = binrule)
#    print(M.facilities[TransList[1]].value)

    M.carbonset = Constraint(expr = summation(M.facilities, M.carbon, index = M.sources) == M.carbonsum)
#    M.Co2limit = Constraint(expr = M.carbonsum <= CO2)  
    
    def objrule(model):
       ob = summation(model.facilities,model.c, index=M.stations) + summation(model.cape, model.isopen, index=M.stations)\
       + (summation(model.locopen, model.loccap, index = model.locations) + summation(model.hydrouse, model.locopex, index = model.locations)) # + summation(model.locopen, model.locopex, index = M.locations) # 
       for i in range(hyn):
           for j in range(locationNum):
               ob = ob + model.hydrouse[CO2LocList[j]]*model.assignments[i*locationNum + j]*costPKGMatrix[i,j]
       return ob

    
    M.Obj = Objective(rule = objrule, sense = minimize)
            
    return M

def opti(model):
    opt = SolverFactory('gurobi', tee = True)
    results = opt.solve(model, tee = True)
    #print(model.display())
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

SourceIn    = pd.read_excel('input.xlsx', 'Sources', index_col=None, na_values=['NA'])
SinkIn      = pd.read_excel('input.xlsx', 'Sinks', index_col=None, na_values=['NA'])
TransIn     = pd.read_excel('input.xlsx', 'Transformers', index_col=None, na_values=['NA'])
HubIn      = pd.read_excel('input.xlsx', 'Hubs', index_col=None, na_values=['NA'])
ConnIn      = pd.read_excel('input.xlsx', 'Connectors', index_col=None, na_values=['NA'])
CO2LocIn  = pd.read_excel('input.xlsx', 'CO2Locations', index_col=None, na_values=['NA'])
RestrIn      = pd.read_excel('input.xlsx', 'Restrictions', index_col=None, na_values=['NA'])

SourceList     = []
SinkList       = []
TransList      = []
H2TransList    = []
HubList        = []
ConnList       = []
CO2LocList     = []
FuelTypeList   = []
DemandTypeList = []
outcols = ['Total Cost', 'CO2']


#Import restrictions
CO2Max = RestrIn.loc[0,'CO2 Max']
wacc = RestrIn.loc[0, 'WACC']
lifetime = RestrIn.loc[0, 'Lifetime']
minCO2PlantSize = RestrIn.loc[0, 'MinCO2PlantSize']
maxCO2PlantSize = RestrIn.loc[0, 'MaxCO2PlantSize']

#Energy sources available from sources
for i in range(len(SourceIn.index)):
    if not SourceIn.loc[i,'EnergyType'] in FuelTypeList:
        FuelTypeList.append(SourceIn.loc[i,'EnergyType'])
        outcols.append(SourceIn.loc[i,'EnergyType'])
        
#Energy types demanded at sinks     
for i in range(len(SinkIn.index)):
    if not SinkIn.loc[i, 'EnergyType'] in DemandTypeList:
        DemandTypeList.append(SinkIn.loc[i, 'EnergyType'])

#All energy types 
EnergyList = FuelTypeList + DemandTypeList

#Initialize the connectors        
for i in range(len(ConnIn.index)):
    ConnList.append(Connection(name = ConnIn.loc[i,'Name'],
                              inp = ConnIn.loc[i,'From'],
                              out = ConnIn.loc[i,'To'],
                              energyType = ConnIn.loc[i,'EnergyType']))

#Initialize the Sources
for i in range(len(SourceIn.index)):
    SourceList.append(Source(name = SourceIn.loc[i,'Name'],
                             energyType = SourceIn.loc[i,'EnergyType'],
                             capex=SourceIn.loc[i,'Capex'], 
                             opex = SourceIn.loc[i,'Opex'], 
                             CO2 = SourceIn.loc[i,'CO2'],
                             minProd = SourceIn.loc[i,'MinProduction'],
                             maxProd = SourceIn.loc[i, 'MaxProduction']))
    
    for con in ConnList:
        if con.inp==SourceList[i].name and con.energyType==SourceList[i].energyType:
            SourceList[i].outcons.append(con)

#Initialize the sinks
for i in range(len(SinkIn.index)):
    SinkList.append(Sink(name = SinkIn.loc[i,'Name'],
                         energyType = SinkIn.loc[i,'EnergyType'],
                         capex = SinkIn.loc[i,'Capex'],
                         opex = SinkIn.loc[i,'Opex'],
                         demand = SinkIn.loc[i,'Demand']))
    
    for con in ConnList:
        if con.out==SinkList[i].name and con.energyType==SinkList[i].energyType:
            SinkList[i].incons.append(con)

#Initialize the transfomers
hyn = 0
for i in range(len(TransIn.index)):
    TransList.append(Transformer(name = TransIn.loc[i,'Name'],
                                 capex = TransIn.loc[i,'Capex'],
                                 opex = TransIn.loc[i,'Opex'],
                                 totalEff = TransIn.loc[i,'TotalEff'],
                                 outMin = TransIn.loc[i, 'OutMin'],
                                 outMax = TransIn.loc[i, 'OutMax'],
                                 CO2Ratio = TransIn.loc[i, 'CO2 Kg/MJ']))
    
    outcols.append(TransList[i].name + 'Production')
    if TransIn.loc[i,'Input0'] == 'hydrogen':
        TransList[i].hynum = hyn
        hyn = hyn + 1
        H2TransList.append(TransList[i])
    
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

#Initialize the Hubs   
for i in range(len(HubIn.index)):
    HubList.append(Hub(name = HubIn.loc[i,'Name'],
                       energyType = HubIn.loc[i,'EnergyType'],
                       capex = HubIn.loc[i,'Capex'],
                       opex = HubIn.loc[i,'Opex']))
    
    outcols.append(HubList[i].name + 'Usage')
    for con in ConnList:
        if con.out==HubList[i].name and con.energyType==HubList[i].energyType:
            HubList[i].incons.append(con)
        elif con.inp==HubList[i].name and con.energyType==HubList[i].energyType:
            HubList[i].outcons.append(con)
    
#Initialize the CO2 Locations
j = 0

for i in range(len(CO2LocIn.index)):                          #Checks to make sure plant capacity is within given bounds
    if(CO2LocIn.loc[i,'Plant size [MW]'] >= minCO2PlantSize): #according to the input restrictions sheet                                                             
        CO2LocList.append(CO2Loc(name = CO2LocIn.loc[i, 'FacilityName'],
                                 ind = j,                               #This if statement only checks mins, maxes
                                 postal = CO2LocIn.loc[i, 'PostalCode'],  #Are checked below, during the calculation
                                 dist = CO2LocIn.loc[i, 'Spalte2'],       #of CO2 location properties by checkMinMax()
                                 cap = CO2LocIn.loc[i, 'Plant size [MW]'],
                                 costPKG = CO2LocIn.loc[i, 'CO2 Cost pkg']))
        j = j + 1  
   
locationNum = j

#Calculate CO2 location properties


for CO2Loc in CO2LocList:
    
    CO2Loc.checkMinMax()
    CO2Loc.findCapex()
    CO2Loc.findDirOpex()
    CO2Loc.findIndOpex()
    CO2Loc.changeCapUnit()
    
    CO2Loc.K = CO2Loc.capex * ((wacc*(wacc + 1)**lifetime)/((wacc+1)**lifetime -1))

checkModel(ConnList, EnergyList)

costPKGMatrix = np.zeros((len(H2TransList),len(CO2LocList)))
for Trans in H2TransList:
    for loc in CO2LocList:
        costPKGMatrix[Trans.hynum,loc.ind] = Trans.CO2Ratio * loc.costPKG

model = createModel(SourceList, SinkList, TransList, ConnList, HubList, CO2LocList, CO2 = CO2Max)

results = opti(model)
#
#for loc in CO2LocList:
#    if model.locopen[loc].value > 10**-12:
#        print(loc.ind)
#        print(model.hydrouse[loc].value)
#        print(model.assignments[loc.ind].value)
#        print(model.assignments[locationNum + loc.ind].value)
#        print(model.locopen[loc].value)

#Output formatting starts here
    
    #Format first sheet

outdf = pd.DataFrame(np.zeros((1,len(outcols))), columns = outcols)
outdf.at[0, 'Total Cost'] = model.Obj()
outdf.at[0, 'CO2'] = model.carbonsum.value

for fac in model.stations:
    if isinstance(fac, Source):
        outdf.at[0, fac.energyType] = model.facilities[fac].value
    elif isinstance(fac, Transformer):
        outdf.at[0, fac.name + 'Production'] = model.facilities[fac].value
        for con in fac.outcons:
            outdf.at[0, fac.name + '-' + con.energyType] = model.connections[con].value
    else:
        outdf.at[0, fac.name + 'Usage'] = model.facilities[fac].value
       
    #Format second sheet

locNames = []
locPosts = []
locAmts = []
locProcs = []

for loc in CO2LocList:
    if model.locopen[loc].value > 0.0000001:
        locNames.append(loc.name)
        locPosts.append(loc.postal)
        locAmts.append(model.hydrouse[loc].value)
        for j in range(hyn):
            if model.assignments[j*locationNum + loc.ind] == 1:
                n = j
                break
            
        procName = H2TransList[n].name
        locProcs.append(procName)

locdf = pd.DataFrame({'Name': locNames,
                      'Postal Code': locPosts,
                      'Amount Used': locAmts,
                      'Process Used': locProcs})    

writer = pd.ExcelWriter('output.xlsx', engine='xlsxwriter')
    
outdf.to_excel(writer, sheet_name='FacilityInfo')
locdf.to_excel(writer, sheet_name='CO2LocationInfo')

writer.save()

checkModel(ConnList, EnergyList)
