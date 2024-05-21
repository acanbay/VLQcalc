"""
Author: Ali Can Canbay - acanbay
Date  : 21 May 2024

This module contains the kappa calculation methods.
calcRatioKappas calculates kappas based on Gamma/Mass ratios using a modified genetic algorithm [1].

====REFERENCES====
[1] Canbay, A. C., (2022). 
    Optimizing KB-mirror Parameters in XRayTracer with Modified Genetic Algorithm. 
    International Hybrid Workshop on "Start-to-End Beamline Optimization for Synchrotron Radiation and 
    Free-Electron Laser Facilities through Artificial Intelligence Approaches", 23-24 November 2022, 
    DESY-Hamburg-Germany.
    https://indico.desy.de/event/36469/

"""
import random, math, cmath
###########################
class optimizeKappas():
    sigma = math.pi/3
    Mass = None
    Kappas = None
    def __init__(self,calcDecay,Mass,BRs,LR):
        self.calcDecay = calcDecay
        self.Mass = Mass
        self.BRs = BRs
        self.LR=LR
        self.nGen = 100
        self.X = self.generatePopulation()
        self.kappas = []
        self.results = []
    
    def function(self,Kappas):
        DH,DW,DZ,Gamma = self.calcDecay(self.Mass,Kappas,self.LR)
        return (DH/Gamma - self.BRs[0])**2 + (DW/Gamma - self.BRs[1])**2 + (DZ/Gamma - self.BRs[2])**2

    def generatePopulation(self):
        X=[]
        for i in range(100):
            chromosome = []
            for val in range(3):
                if self.BRs[val] == 0:
                    chromosome.append(0)
                else:
                    chromosome.append( random.uniform(0,1) )
            X.append(chromosome[:])
        return X
    
    def rouletteChoice(self, pool):
        index = abs( int( random.gauss(0,10) ) )
        return pool[index]

    def crossover(self, parent1, parent2):
        offSpring = []
        for i in range(3):
            if self.BRs[i] == 0:
                offSpring.append(0)
            else:
                if random.random() < 0.33:
                    offSpring.append( (parent1[i]+parent2[i])/2 )
                else:
                    if random.random() < 0.5:
                        offSpring.append( parent1[i] )
                    else:
                        offSpring.append( parent2[i] )
        return offSpring
    
    def mutation(self, offSpring):
        for i in range(3):
            if self.BRs[i] == 0:
                offSpring.append(0)
            else:
                if random.random() < 0.33:
                    gene = offSpring[i] * random.gauss(1, self.sigma)
                    if  0 < gene < 1:
                        offSpring[i] = gene
        return offSpring

    def nextGeneration(self):
        population = self.X[:]
        population=sorted(population, key=self.function)

        new_population=population[:1]
        for i in range(99):
            parent1 = self.rouletteChoice(population)
            parent2 = self.rouletteChoice(population)
            offSpring = self.crossover( parent1, parent2 )
            offSpring = self.mutation(offSpring)
            new_population.append(offSpring)
        self.X=new_population

        if len(self.results)>1:
            dx = math.sqrt( sum( [ (self.kappas[-2][i] - self.kappas[-1][i])**2 for i in range(3) ] ) )
            dy = self.results[-2] - self.results[-1] 
            if dx != 0 and dy>0:
                deriv = dy/dx
                self.sigma=deriv/3
        
        if len(self.results)>1:
            dx = math.sqrt( sum( [ (self.kappas[-2][i] - self.kappas[-1][i])**2 for i in range(3) ] ) )
            dy = self.results[-2] - self.results[-1] 
            if dx != 0 and dy>0:
                deriv = dy/dx
                if deriv > math.pi/3:
                    deriv = math.pi/3
                self.sigma=deriv

    def optimize(self):
        for n in range(self.nGen):
            self.nextGeneration()
            self.kappas.append( self.X[0] )
            self.results.append( self.function( self.kappas[-1] ) )
        return self.kappas[-1]

def calcRatioKappas(VLQ):
    Kappas=[]
    if VLQ.type in ["VLX","VLY"]:
        for Mass in VLQ.Mass:
            Gamma = VLQ.calcDecay(Mass, 0.5, VLQ.LR)
            Kappas.append( round( cmath.sqrt( 0.5**2 * VLQ.Ratio / ( Gamma / Mass ) ).real, 6 ))
        VLQ.KappaW = Kappas
    else:
        Kappas=[[],[],[]]
        for Mass in VLQ.Mass:
            ga = optimizeKappas(VLQ.calcDecay,Mass,VLQ.BRs,VLQ.LR)
            ks = ga.optimize()
            ## Normalizing
            sumks = sum(ks)
            ks = [k/sumks for k in ks ]
            Kappas[0].append(ks[0])
            Kappas[1].append(ks[1])
            Kappas[2].append(ks[2])
        for i in range(len(VLQ.Mass)):
            DH, DW, DZ, Gamma = VLQ.calcDecay( VLQ.Mass[i], [Kappas[0][i], Kappas[1][i], Kappas[2][i]],VLQ.LR)
            Kappas[0][i] = round( (cmath.sqrt( Kappas[0][i]**2 * VLQ.Ratio / ( Gamma / VLQ.Mass[i] ) ).real), 6)
            Kappas[1][i] = round( (cmath.sqrt( Kappas[1][i]**2 * VLQ.Ratio / ( Gamma / VLQ.Mass[i] ) ).real), 6)
            Kappas[2][i] = round( (cmath.sqrt( Kappas[2][i]**2 * VLQ.Ratio / ( Gamma / VLQ.Mass[i] ) ).real), 6)
        VLQ.KappaH = Kappas[0]
        VLQ.KappaW = Kappas[1]
        VLQ.KappaZ = Kappas[2]