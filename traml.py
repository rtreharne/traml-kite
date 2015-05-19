from tabulate import tabulate
from numpy import *
import os
from glob import glob
from csv import reader
from pylab import *

class Stack:


    def __init__(self, filename=None):
        self.BASE = os.path.dirname(os.path.abspath(__file__))+'/library/'
        self.config = []
        self.substrate('OW') 
        self.create_library()

    def plot(self, data, twin='off'):
        x = list(zip(*data)[0]) 
        y = list(zip(*data)[1])
	plot(x, y)
	show()
    
    def film_data(self, film_id):
        filename = self.library_dict[self.config[film_id][0]]
        with open(filename, 'rU') as f:
	    csvdata = reader(f)
	    next(csvdata) # Skip first row
	    rows = list(csvdata)
	    float_rows = [[float(i) for i in row] for row in rows] # Convert all data to float
        return float_rows 

    def create_library(self):
        file_list = sort(os.listdir(self.BASE))
        dict_list = []
        for i in range(len(file_list)):
            dict_list.append([file_list[i].replace('.csv', ''), self.BASE+file_list[i]])
        self.library_dict = dict(dict_list)
 
    def library(self):
        for key, value in self.library_dict.iteritems():
            print key

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
	elif len(self.config)>1:
	   del self.config[loc]
	
	self.table()

    def table(self):
        table = []
        for i in range(len(self.config)):
            table.append(self.config[i][:])
            table[i].insert(0, i)

        print tabulate(table,
	               headers=['#', 'Material', 'Thickness (nm)', 'Type'], 
		       tablefmt='orgtbl')






