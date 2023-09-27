#!/usr/bin/env python
# coding: utf-8

# In[1]:


import abc
import numpy as np
from enum import Enum, IntEnum
from datetime import datetime
import time
import math
import datetime
import matplotlib.pyplot as plt
import ctypes
from tkinter import *
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)


# # Global Function(Identify the date)

# In[ ]:


def days_between(d1,d2):
    '''
    d1 = start date
    d2 = end date
    This function provided differences between dates
    '''
    d1 = datetime.datetime.strptime(d1, "%d/%m/%Y") #start date
    d2 = datetime.datetime.strptime(d2, "%d/%m/%Y") #end date
    return (d2 - d1).days

def days_after(d1,d2):
    '''
    This equation determine whether the relation between dates
    '''
    newdate1 = time.strptime(d1, "%d/%m/%Y")
    newdate2 = time.strptime(d2, "%d/%m/%Y")
    determine = newdate1 <= newdate2 
    if determine == True:
        response_after = 1    #d2 is after d1
    else:
        response_after = 0    #d2 is before d1
    return response_after

def days_before(d1,d2):
    '''
    This equation determine whether the relation between dates
    '''
    newdate1 = time.strptime(d1, "%d/%m/%Y")
    newdate2 = time.strptime(d2, "%d/%m/%Y")
    determine = newdate1 >= newdate2 
    if determine == True:
        response_before = 1    #d2 is before d1
    else:
        response_before = 0    #d2 is after d1
    return response_before


# # Interpolation Function for lambda

# In[185]:


class Interpolation():
    def __init__(self, start_date, end_date, today, ds_rate):
        self.start = start_date
        self.end = end_date
        self.quote = quote_value
        self.today = today
        self.length = len(start_date)
        self.ds_rate = ds_rate
        
        if (not self.valid()):
            raise Exception("Size of start date and end date should be the same")
        
    def getLength(self):
        return self.length
    
    def valid(self):
        if (len(self.start)!=len(self.end)):
            return False
        return True
    
    def z_multi(self):
        '''
        Z quote for multi instruments
        '''
        ds_rate_n = len(self.ds_rate) #assume it has len = 10 (0-9)
        quote = []             #0-10
        quote += [0]
        index = list(range(1,ds_rate_n+1))            #1-10
        for i in index:
            quote += [self.ds_rate[i-1]*0.01]
        return quote #it should be a n+1  
    
    def solve_single_Q_test(self, t): #aim to create a 1xn Q matrix
        '''
        This function is only designed for solving the single t which will only be used in the final cap factor
        after calculating the value of lambda
        '''
        n_col = len(self.start)
        matrix_Q_single = np.zeros(n_col)
        for ii in range(n_col): #1-n column
            response_after = days_after(self.start[ii],t) #1 for true, 0 for false
            response_before = days_before(self.end[ii],t) #1 for true, 0 for false
            if response_after == 1 and response_before == 1:
                diff = days_between(self.start[ii], t)
                matrix_Q_single[ii] = 1/3*(diff)**3
                        
            elif response_after == 0: #before the si
                matrix_Q_single[ii] = 0
                        
            elif response_before == 0: #after the ei
                diff_start_end = days_between(self.start[ii], self.end[ii])
                diff_end_t = days_between(self.end[ii], t)
                matrix_Q_single[ii] = 1/3 * (diff_start_end**3) +                                    (diff_start_end**2)*diff_end_t +                                    diff_start_end * (diff_end_t**2)
        return matrix_Q_single  #1xn matrix of Q
    
      
    def solve_Q(self): #aim to create a nxn Q matrix
        '''
        Create the matrix of Q with size nxn, the matrix contain the Q for t_1...t_n and with start_date_1 to n
        and end_date_1 to n
        '''
        if (not self.valid()):
            raise Exception("Size of start dates and end dates should be the same")
        
        else:
            n = len(self.start)
            matrix_Q = np.zeros((n,n))
            for i in range(n): #1-n row, also represent ti
                for ii in range(n): #1-n column
                    response_after = days_after(self.start[ii],self.end[i]) #1 for true, 0 for false
                    response_before = days_before(self.end[ii],self.end[i]) #1 for true, 0 for false
                    if response_after == 1 and response_before == 1:
                        diff = days_between(self.start[ii], self.end[i])
                        matrix_Q[i][ii] = 1/3 * (diff)**3
                        
                    elif response_after == 0: #before the si
                        matrix_Q[i][ii] = 0
                        
                    elif response_before == 0: #after the ei
                        diff_start_end = days_between(self.start[ii], self.end[ii])
                        diff_start_xValues = days_between(self.end[ii], self.end[i])
                        matrix_Q[i][ii] = 1/3 * (diff_start_end**3) +                                          (diff_start_end**2)*diff_start_xValues +                                          diff_start_end * (diff_start_xValues**2)
        return matrix_Q  #nxn matrix of Q
        
    def matrix(self):
        '''
        Create the matrix with size (n+1)x(n+1), and the matrix is eventually used to solve lambda
        '''
        Q_t = self.solve_Q() #Q_t should be a matrix containing Q1,Q2,...,Qn >>>> (nxn matrix)
        #print(Q_t)
        n = len(self.end)
        matrix = [] #create a (n+1)x(n+1) empty matrix
        first_row = []

        for es in range(n+1):
            if es == 0:
                first_row += [0]
            else:
                diff = days_between(self.start[es-1], self.end[es-1])
                diff = -diff*2 #test
                first_row += [diff] #insert first row
        matrix += [first_row]
        #print(matrix)
        
        for row in range(n): #row 1-n
            rows = []
            for col in range(n+1): #column 1-(n+1)
                if col == 0:
                    rows += [1]
                else:
                    t_i = days_between(self.today, self.end[row])
                    rows += [Q_t[row][col-1]/t_i] 
            matrix += [rows]
                
        return matrix #should be a (n+1)x(n+1) matrix
    
    def solve_lambda_multi(self):
        '''
        Calculate the value of lambda from lambda0,lambda1,...,lambda_n
        Lambda for multi instruments 
        '''
        matrix = self.matrix()
        quote = self.z_multi()
        #lambda_list = np.linalg.inv(matrix)*quote
        invert = np.linalg.inv(matrix)
        #print(invert)
        lambda_list = np.dot(np.linalg.inv(matrix),quote)
        return lambda_list


# # Quadratic Forward Function

# In[ ]:


class QuadraticInterpolator(Interpolation):
    def __init__(self, start_date, end_date, today, ds_rate):
        Interpolation.__init__(self, start_date, end_date, today, ds_rate)    
    
    def solve_zc_multi(self,t):
        '''
        zc for multi instruments
        xValue = list of unknown time we used to solve for lambda
        t = the unknown time we want to solve using interpolation
        quote_value = market quote, corresponding to the xValue
        '''
        #determine whether t is interpolate or extrapolate
        determine_inter = days_before(self.start[0],t) #if before the first
        determine_inter2 = days_after(self.end[-1],t) #if after the last
        matrix = self.matrix()
        lambda_value = self.solve_lambda_multi()
        Q_t_single = self.solve_single_Q_test(t)
        t_ori = days_between(self.today, t)
        #print("the lambda value is", lambda_value)
        #print("the Q is", Q_t_single)
        N = len(lambda_value)
        if determine_inter == 1: #t is before the first
            zc = 0
        elif determine_inter2 == 1: #t is after the last
            zc = lambda_value[0]
            value = 0
            for i in range(N-1):
                add = (Q_t_single[i]*lambda_value[i+1])/t_ori
                value = value + add
            zc = value + zc
        else: #t is in the period
            zc = lambda_value[0]
            #print("the ln factor is", log_cap_factor)
            value2 = 0
            for ii in range(N-1):
                add2 = (Q_t_single[ii]*lambda_value[ii+1])/t_ori
                value2 = value2 + add2
            zc = zc + value2
        df = math.exp(-zc*t_ori/365)
        zc_final = zc*100
        
        number = 0
        window = Tk()
        window.title("Programme")
        window.geometry('350x250')

        label = Label(window, text=number)
        label.grid(column=0,row=0)

        #ttt=zc*100
        ttt = df
        label.config(text=ttt)

        button = Button(window,text="Value of discount factor for current date")
        button.grid(column=1, row=2)

        window.mainloop()
            
        #return df      
        #return zc*100
        
    def solve_zc_multi_data(self,t):
        '''
        zc for multi instruments
        xValue = list of unknown time we used to solve for lambda
        t = the unknown time we want to solve using interpolation
        quote_value = market quote, corresponding to the xValue
        '''
        #determine whether t is interpolate or extrapolate
        determine_inter = days_before(self.start[0],t) #if before the first
        determine_inter2 = days_after(self.end[-1],t) #if after the last
        matrix = self.matrix()
        lambda_value = self.solve_lambda_multi()
        Q_t_single = self.solve_single_Q_test(t)
        t_ori = days_between(self.today, t)
        #print("the lambda value is", lambda_value)
        #print("the Q is", Q_t_single)
        N = len(lambda_value)
        if determine_inter == 1: #t is before the first
            zc = 0
        elif determine_inter2 == 1: #t is after the last
            zc = lambda_value[0]
            value = 0
            for i in range(N-1):
                add = (Q_t_single[i]*lambda_value[i+1])/t_ori
                value = value + add
            zc = value + zc
        else: #t is in the period
            zc = lambda_value[0]
            #print("the ln factor is", log_cap_factor)
            value2 = 0
            for ii in range(N-1):
                add2 = (Q_t_single[ii]*lambda_value[ii+1])/t_ori
                value2 = value2 + add2
            zc = zc + value2
        df = math.exp(-zc*t_ori/365)
        #print(df)
        zc_final = zc*100
        
        return zc_final  
    
    def check_with_mx(self, t_series):
        n = len(t_series)
        check_t = []
        for i in range(n):
            zc_t = self.solve_zc_multi(t_series[i])
            check_t += [zc_t]
            
        return check_t
    
    def tk_zc_plot(self, nod_t):
        window = Tk()
        window.title('Plotting in Tkinter')
        window.geometry("500x500")
        fig = Figure(figsize = (5, 5), dpi = 100)
        zc_list = []
        len_t = len(nod_t)
        for i in range(len_t):
            nod= self.solve_zc_multi_data(nod_t[i])
            zc_list += [nod]
        x_axis = range(1,len_t)
        plot1 = fig.add_subplot(111)
        plot1.plot(x_axis,zc_list[1:len_t])
        canvas = FigureCanvasTkAgg(fig, master = window)  
        canvas.draw()
        canvas.get_tk_widget().pack()
        toolbar = NavigationToolbar2Tk(canvas,window)
        toolbar.update()
        canvas.get_tk_widget().pack()
        plot_button = Button(master = window, 
                     height = 2, 
                     width = 10,
                     text = "Plot")
        #plot_button.pack()
        window.mainloop()
        
    def check_sensitivity(self,t,flows):
        matrix = self.matrix()
        #print("the matrix is",np.linalg.inv(matrix))
        n = len(matrix)
        Q = self.solve_single_Q_test(t)
        #print("the Q is",Q)
        Q_t = []
        for i in range(n):
            if i == 0:
                Q_t += [1]
            else:
                t_i = days_between(self.today, t)
                Q_t += [Q[i-1]/t_i]     
        test_sensitivity = np.dot(Q_t,np.linalg.inv(matrix))
        time_period = days_between(self.today,t)
        z = self.solve_zc_multi_data(t)/100
        flow_sum = sum(flows)
        sensitivity = (test_sensitivity*math.exp(-z*time_period/365)*-time_period/365*(flow_sum))/10000

        return sensitivity[1:]


# # Test database

# In[ ]:


start_date_2 = ["07/12/2022","07/12/2022","07/12/2022","20/03/2024","18/06/2025","17/12/2025"]
end_date_2 = ["07/01/2023","07/06/2023","07/12/2023","20/06/2024","18/09/2025","17/03/2026"]
ds_rate = [0.999999999999896,1.09999999999994,1.20000000000001,2.03661722517255,2.79066733721628,2.8406346644081]
#quote_value_2 = [1,1.1,1.2,4.265,3.16,3.08499999999999]
today = "07/12/2022"
t = "05/10/2023"
t2 = "23/04/2017"
t3 = "23/04/2015"
inst_1 = "09/06/2023"
i = "07/06/2023"
flows = [-5053823.5925004,-1000000000]


# In[ ]:


'''
create time period from 07/12/2022 to 07/12/2026
'''
start = datetime.datetime.strptime("07/12/2022", "%d/%m/%Y")
end = datetime.datetime.strptime("07/12/2026", "%d/%m/%Y")
date_array = (start + datetime.timedelta(days=x) for x in range(0, (end-start).days))
time_t = []
for date_object in date_array:
    time_t += [date_object.strftime("%d/%m/%Y")]


# In[ ]:


curveB = QuadraticInterpolator(start_date_2, end_date_2, quote_value_2, today, ds_rate)
curveB.solve_lambda_multi()


# In[ ]:


curveB.check_sensitivity(inst_1,flows)


# In[ ]:


curveB.solve_zc_multi(i)


# In[ ]:


curveB.tk_zc_plot(time_t)


# In[ ]:


curveB.matrix()


# In[ ]:


curveB.solve_lambda_multi()


# In[ ]:


y = curveB.check_with_mx(xValues_2, time_t)


# In[ ]:




