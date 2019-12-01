# -*- coding: utf-8 -*-
"""
Created on Thr Nov 29 2019

@author: bagheri

using Genetic Algorithm for finding minimum of f(x,y)
"""
import numpy as np
import random
import matplotlib.pyplot as plt
from statistics import mean

class GeneticAlgorithm():
    def __init__(self,problem):
        self.problem=problem
        self.generation_fitness=[]
        self.generation=[]
        self.history=[]
        
    def findOpt(self,genomSize,min_val,max_val,max_iters,gen_pop,cross_prob,mut_prob,precision,error,maximize):
        self.min_val=min_val
        self.max_val=max_val
        self.gen_pop=gen_pop
        self.cross_prob=cross_prob
        self.genomSize=genomSize
        self.mut_prob=mut_prob
        self.precision=precision
        self.maximize=maximize
        
        self.generation.clear()
        self.history.clear()
        self.oldfitness=[]
      
        self.initializePopulation()
        
        converge=False
        itr=0

        while itr<max_iters and not(converge) :
            self.oldfitness.clear()
            self.oldfitness=self.generation_fitness.copy()
            selectedGenoms=self.selection()
            self.crossOver(selectedGenoms)
            self.mutation()
            self.calGenerationFitness()
            itr=itr+1
            if(self.maximize):
                if abs(self.history[-1][1]-self.history[-2][1])<error: converge=True
                else: converge=False
            else:
                if abs(self.history[-1][0]-self.history[-2][0])<error: converge=True
                else: converge=False
 
        
        return self.history,itr
      
    def initializePopulation(self):
        for i in range(0,self.gen_pop):
            g=[]
            g.clear()
            for j in range(0,self.genomSize):
                x=format(random.randint(0,pow(2,self.precision)-1),'b')
                if len(x)<self.precision:
                    x=(self.precision-len(x))*'0'+x
                g.append(x)   
            self.generation.append(g)
       
        self.calGenerationFitness()
        
    def calGenerationFitness(self):
        self.generation_fitness.clear()
        for genom in self.generation:
            f=self.calGenomFitness(genom)
            self.generation_fitness.append(f)
        self.history.append([min(self.generation_fitness),max(self.generation_fitness),mean(self.generation_fitness)])
    

    def calGenomFitness(self,genom):
        state=[]
        for gen in genom:
            x=int(gen,2)/(pow(2,self.precision)-1)*(self.max_val-self.min_val)
            state.append(x)
        f=self.problem(state)
        return f
    
    def selection(self):
        selectedGenoms=[]
        self.generation_fitness.sort(reverse=self.maximize)
        for i in range(0,int(abs(self.gen_pop*self.cross_prob))):
            for s in self.generation:
                
                if self.generation_fitness[i]==self.calGenomFitness(s): 
                    selectedGenoms.append(s)
                    break
        return selectedGenoms
                

    def crossOver(self,parents):
        no_offSprings=0
        newGeneration=[]
        tempGen=[]
        while no_offSprings<self.gen_pop:
            if len(parents)==0:
                parents=tempGen.copy()
                tempGen.clear()
            xx=[x for x in parents if x not in tempGen]
            if len(xx)==0:
                parents=tempGen.copy()
                tempGen.clear()
                xx=[x for x in parents if x not in tempGen]
            p1=random.choice(xx)
            p2=random.choice(xx)
            while (p1==p2) :
                p2=random.choice(parents)
            tempGen.append(p1)
            tempGen.append(p2)
            parents.remove(p1)
            parents.remove(p2)
            fraction=random.randint(0,self.precision-1)
            if fraction==0: fraction=1
            c1=[]
            c2=[]
            for g in range(0,self.genomSize):
               c1.append(p1[g][0:fraction]+p2[g][fraction:])
               c2.append(p2[g][0:fraction]+p1[g][fraction:])
            
            newGeneration.append(c1)
            newGeneration.append(c2)
            no_offSprings=no_offSprings+2
        self.generation.clear()
        self.generation=newGeneration

  
    def mutation(self):
        for genom in self.generation:
            prob=random.random()
            for gen in genom:
                if prob<self.mut_prob:
                    i=int(random.uniform(0,self.precision))
                    if gen[i]=='1': 
                        gen=gen[:i]+'0'+gen[i+1:]
                    else: 
                        gen=gen[:i]+'1'+gen[i+1:]
    
    
    
def myProblem(state):
    return ((state[0]-3)*(state[0]-3))+((state[1]-2)*(state[1]-2))


def main():
    ga=GeneticAlgorithm(myProblem)
    history,itr=ga.findOpt(genomSize=2,min_val=0,max_val=6,max_iters=10,gen_pop=10,cross_prob=0.8,mut_prob=.05,precision=8,error=0.01,maximize=False)
    minf=[]
    maxf=[]
    avef=[]
    for h in history:
        minf.append(h[0])
        maxf.append(h[1])
        avef.append(h[2])
    
    plt.plot(range(len(minf)), minf, label='min')
    plt.plot(range(len(maxf)), maxf, label='max')
    plt.plot(range(len(avef)), avef ,label='average')
    plt.title('')
    plt.xlabel('Iteration')
    plt.ylabel('Fittnes')
    plt.legend()
    plt.show()

main()


