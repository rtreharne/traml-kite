from tabulate import tabulate
from numpy import *
from copy import copy, deepcopy

class Stack:


    def __init__(self, filename=None):
        self.substrate('glass') 
        self.config = []

    def substrate(self, name, d='--', film_type='substrate'):
        self.config.append([name, d,  film_type])

    def add(self, name, d=100, film_type='passive', loc=None):
        if loc==None:
	    self.config.append([name, d, film_type])
	else:
	    if loc>0:
	        self.config.insert(loc, [name, d, film_type])

	self.table()

    def remove(self, loc=None):
        if loc==None and len(self.config)>1:
	   self.config.pop()
	else:
	   del self.config[loc]
	
	self.table()

    def table(self):
        table = []
        for i in range(len(self.config)):
            table.append(self.config[i-1][:])
            table[i-1].insert(0, i)

        print tabulate(table,
	               headers=['#', 'Material', 'Thickness (nm)', 'Type'], 
		       tablefmt='orgtbl')





