#!/usr/bin/python
# coding=utf-8

#Created by Parth Patel, DBI @ University of Delaware, Newark, Delaware 19717
#Date created: 11/23/2015, Updated:08/14/18

## This script calculates the usage  position-specific base usage of  between two types of smallRNAs (i.e., phasiRNAs and  non-phasiRNAs(P4-siRNAs))

import os
import numpy as np
import scipy.stats as sp
import matplotlib as mpl
import sys
## agg backend is used to create plot as a .png file
mpl.use('agg')
#mpl.use('TkAgg') for displaying plt.show() to show figures on windows
from pylab import *
import matplotlib.pyplot as plt
#import beanplot
import math
from scipy import stats
from scipy.stats import mannwhitneyu
from scipy.stats import ranksums
import subprocess

# (*) To communicate with Plotly's server, sign in with credentials file
from plotly import tools
import plotly.plotly as py
# (*) Graph objects to piece together plots
from plotly.graph_objs import *
import plotly.graph_objs as go
from plotly.graph_objs import Scatter, Data, Layout
import collections
import plotly



########################## Sign in to plotly  ###############################


py.sign_in("ParthPatel", "56xws1snxp")
plotly.tools.set_credentials_file(username='ParthPatel', api_key='56xws1snxp')



########################## GLOBAL VARIABLES ###############################

positiveSet= sys.argv[1] #"phasiRNAs.txt" 
negativeSet=sys.argv[2] #"hc-siRNAs.txt"
size=24 # size of sRNAs in both files
p_val=0.001 # p_value to determine the signficance of the result


##############################################################################################################
################################# Draw position specific base usage###########################################
##############################################################################################################


def drawPFM(seq_Comp,seq_comp_non_phasi,n,significant_positions):
    
    ## a stacked bar plot 
    
    py.sign_in("ParthPatel", "56xws1snxp")
    
    x=np.linspace(1,n,n) # For example, create number between 1 to 21 21 times (meaning with step of 1)
    y1= seq_Comp [0][:] # A  
    y2=seq_Comp [1][:] # C
    y3=seq_Comp [2][:] # G
    y4=seq_Comp [3][:]# U
    
    
    print significant_positions
    
    A_x=[]
    A_y=[]
    C_x=[]
    C_y=[]
    G_x=[]
    G_y=[]
    U_x=[]
    U_y=[]
    
    for i in range(0,len(significant_positions)):
        letter=significant_positions[i][0]
        position=int(significant_positions[i][1:])-1
#        print position,letter
        if letter=='A':
            A_x.append(significant_positions[i][1:])
            A_y.append(y1[position])
        elif letter=='C':
            C_x.append(significant_positions[i][1:])
            C_y.append(y2[position])
        elif letter=='G':
            G_x.append(significant_positions[i][1:])
            G_y.append(y3[position])
        elif letter=='T':
            U_x.append(significant_positions[i][1:])
            U_y.append(y4[position])
            
    
    

    trace1 = go.Scatter(
    x = x,
    y = y1,
    mode = 'lines+markers',
    line=Line(
        color='rgb(199,21,133)',
        dash='solid',
        width=3
    ),
    marker=Marker(color='rgb(102, 0, 51))',symbol='circle-open-dot'),
    name = 'A'
    )
    

    
    trace2 = go.Scatter(
    x = x,
    y = y2,
    mode = 'lines+markers',
    line=Line(
        color='rgb(32,178,170)',
        dash='solid',
        width=3
    ),
    marker=Marker(color='rgb(0, 51, 102)',symbol='circle-open-dot'),
    name = 'C'
    )
    
    trace3 = go.Scatter(
    x = x,
    y = y3,
    mode = 'lines+markers',
     line=Line(
        color='rgb(255, 175, 0)',
        dash='solid',
        width=3
    ),
    marker=Marker(color='rgb(255, 128, 0)',symbol='circle-open-dot'),
    name = 'G'
    )
    
    trace4 = go.Scatter(
    x = x,
    y = y4,
    mode = 'lines+markers',
    line=Line(
        color='rgb(0,186,56)',
        dash='solid',
        width=3
    ),
    marker=Marker(color='rgb(0,102,0)',symbol='circle-open-dot'),
    name = 'U'
    )

    data=[trace1,trace2,trace3,trace4]
    
    layout = Layout(
        title='',
          xaxis=XAxis(
          title='Position',
          ticks='outside',
          tickcolor='rgb(102, 102, 102)',
            tickfont=Font(
                size=14,
                color='rgb(107, 107, 107)'
            )
        ),
        yaxis=YAxis(
            title='Frequency',
            ticks='outside',
            tickcolor='rgb(102, 102, 102)',
            tick0=0,
            dtick=.1,
            range=[0,0.7],
            titlefont=Font(
                size=16,
                color='rgb(107, 107, 107)'
            ),
            tickfont=Font(
                size=14,
                color='rgb(107, 107, 107)'
            )
        ),

        plot_bgcolor='#ffffff',
    )    
      
    annotations=[]

    for i in range (0,len(A_x)):
        annotations.append(dict(x=A_x[i],y=A_y[i],xref='x',yref='y',showarrow=False,text="",arrowhead=7,bordercolor='rgb(102, 0, 51)',bgcolor='rgb(102, 0, 51)',borderwidth=2.3))
    for i in range (0,len(C_x)):
        annotations.append(dict(x=C_x[i],y=C_y[i],xref='x',yref='y',showarrow=False,text="",arrowhead=7,bordercolor='rgb(0, 51, 102)',bgcolor='rgb(0, 51, 102)',borderwidth=2.3))
    for i in range (0,len(G_x)):
        annotations.append(dict(x=G_x[i],y=G_y[i],xref='x',yref='y',showarrow=False,text="",arrowhead=7,bordercolor='rgb(255, 128, 0)',bgcolor='rgb(255, 128, 0)',borderwidth=2.3))
    for i in range (0,len(U_x)):
        annotations.append(dict(x=U_x[i],y=U_y[i],xref='x',yref='y',showarrow=False,text="",arrowhead=7,bordercolor='rgb(0,102,0)',bgcolor='rgb(0,102,0)',borderwidth=2.3))
           
    layout['annotations'] = annotations
    
   
    
    y5= seq_comp_non_phasi [0][:] # A    
    y6=seq_comp_non_phasi [1][:] # C
    y7=seq_comp_non_phasi [2][:] # G
    y8=seq_comp_non_phasi [3][:]# U
    
    # re-initialize it
    A_x=[]
    A_y=[]
    C_x=[]
    C_y=[]
    G_x=[]
    G_y=[]
    U_x=[]
    U_y=[]
    
    for i in range(0,len(significant_positions)):
        letter=significant_positions[i][0]
        position=int(significant_positions[i][1:])-1
        if letter=='A':
            A_x.append(significant_positions[i][1:])
            A_y.append(y5[position])
#            print "A",significant_positions[i][1:],y5[position]
        elif letter=='C':
            C_x.append(significant_positions[i][1:])
            C_y.append(y6[position])
#            print "C",significant_positions[i][1:],y6[position]
        elif letter=='G':
            G_x.append(significant_positions[i][1:])
            G_y.append(y7[position])
#            print "G",significant_positions[i][1:],y7[position]
        elif letter=='T':
            U_x.append(significant_positions[i][1:])
            U_y.append(y8[position])
#            print "T",significant_positions[i][1:],y8[position]
    
    
    
    trace5 = go.Scatter(
    x = x,
    y = y5,
    mode = 'lines+markers',
    line=Line(
        color='rgb(199,21,133)',
        dash='solid',
        width=3
    ),
    marker=Marker(color='rgb(102, 0, 51))',symbol='circle-open-dot'),
    name = 'A'
    )
    
    trace6 = go.Scatter(
    x = x,
    y = y6,
    mode = 'lines+markers',
    line=Line(
        color='rgb(32,178,170)',
        dash='solid',
        width=3
    ),
    marker=Marker(color='rgb(0, 51, 102)',symbol='circle-open-dot'),
    name = 'C'
    )
    
    trace7 = go.Scatter(
    x = x,
    y = y7,
    mode = 'lines+markers',
     line=Line(
        color='rgb(255, 175, 0)',
        dash='solid',
        width=3
    ),
    marker=Marker(color='rgb(255, 128, 0)',symbol='circle-open-dot'),
    name = 'G'
    )
    
    trace8 = go.Scatter(
    x = x,
    y = y8,
    mode = 'lines+markers',
    line=Line(
        color='rgb(0,186,56)',
        dash='solid',
        width=3
    ),
    marker=Marker(color='rgb(0,102,0)',symbol='circle-open-dot'),
    name = 'U'
    )
    
    data1=[trace5,trace6,trace7,trace8]
      
    
    
    layout1 = Layout(
        title= '',
          xaxis=XAxis(
          title='Position',
          ticks='outside',
          tickcolor='rgb(102, 102, 102)',
            tickfont=Font(
                size=14,
                color='rgb(107, 107, 107)'
            )
        ),
        yaxis=YAxis(
            title='Frequency',
            ticks='outside',
            tickcolor='rgb(102, 102, 102)',
            tick0=0,
            dtick=.1,
            range=[0,0.7],
            titlefont=Font(
                size=16,
                color='rgb(107, 107, 107)'
            ),
            tickfont=Font(
                size=14,
                color='rgb(107, 107, 107)'
            )
        ),

        plot_bgcolor='#ffffff',

    )
    
    
    annotations=[]

    for i in range (0,len(A_x)):
        annotations.append(dict(x=A_x[i],y=A_y[i],xref='x',yref='y',showarrow=False,text="",arrowhead=7,bordercolor='rgb(102, 0, 51)',bgcolor='rgb(102, 0, 51)',borderwidth=2.3))
    for i in range (0,len(C_x)):
        annotations.append(dict(x=C_x[i],y=C_y[i],xref='x',yref='y',showarrow=False,text="",arrowhead=7,bordercolor='rgb(0, 51, 102)',bgcolor='rgb(0, 51, 102)',borderwidth=2.3))
    for i in range (0,len(G_x)):
        annotations.append(dict(x=G_x[i],y=G_y[i],xref='x',yref='y',showarrow=False,text="",arrowhead=7,bordercolor='rgb(255, 128, 0)',bgcolor='rgb(255, 128, 0)',borderwidth=2.3))
    for i in range (0,len(U_x)):
        annotations.append(dict(x=U_x[i],y=U_y[i],xref='x',yref='y',showarrow=False,text="",arrowhead=7,bordercolor='rgb(0,102,0)',bgcolor='rgb(0,102,0)',borderwidth=2.3))
           
    layout1['annotations'] = annotations
       
    fig = Figure(data=data, layout=layout)
    py.image.save_as(fig,'phasiRNAs',format='png',scale=6)
        
    fig = Figure(data=data1, layout=layout1)
    py.image.save_as(fig,'hc-siRNAs',format='png',scale=6)
    
 


def RankSumTest(group1,group2):

    # Compute the Wilcoxon rank-sum statistic for two samples x and y.
    u, p_value = ranksums(group1, group2)

    return p_value
   

def main(phasiRNA_sequences,non_phasiRNA_sequences,k_mer_list,n,significance):
    
   
   
##################################################################################
##
##   Initialization of all variables
##   This routine counts freuncies of  A/C/G/T at each position
##
##################################################################################   
    
    # num_lines- get umber of lines aka. number of sequences in the files.    
    num_lines = sum(1 for line in open(phasiRNA_sequences))
    num_lines1 = sum(1 for line in open(non_phasiRNA_sequences))
   
    #for phasiRNAs do the below
    
    
    seq_Comp_phasiRNA=np.zeros((4,n),int) #a 4x21 array of int zeros    
    
    #read files
    input=open(phasiRNA_sequences, 'r')
    
  
    for lines in input:
        
        #print "Reading Line: %s" % (lines)
        sequence = lines.strip()
        if not lines: break

        for y in range (0,n):
            if (sequence [y] == 'A'):
                seq_Comp_phasiRNA [0][y] =  seq_Comp_phasiRNA [0][y]+1
            if (sequence [y] == 'C'):
                seq_Comp_phasiRNA [1][y] = seq_Comp_phasiRNA [1][y]+1
            if (sequence [y] == 'G'):
                seq_Comp_phasiRNA [2][y] = seq_Comp_phasiRNA [2][y]+1
            if (sequence [y] == 'T'):
                seq_Comp_phasiRNA [3][y] = seq_Comp_phasiRNA [3][y]+1   
       
 
########################################################################################
##
##   Initialization of all variables
##   This routine gets normalized frequency of nucleotides at each position in phasiRNAs.
##
########################################################################################    
           
    (rows,columns)=seq_Comp_phasiRNA.shape 
    #print  seq_Comp_phasiRNA.shape          
    Norm_Count_phasiRNA=np.zeros((4,columns),float) #a 4x21 array of int zeros
    
    for x in range(rows):
        for y in range(columns):
          
            Norm_Count_phasiRNA[0][y] = float(seq_Comp_phasiRNA[0][y])/num_lines #Total_A's in first position/Total_phasiRNA
            Norm_Count_phasiRNA[1][y] = float(seq_Comp_phasiRNA[1][y])/num_lines #Total_C's in first position/Total_phasiRNA
            Norm_Count_phasiRNA[2][y] = float(seq_Comp_phasiRNA[2][y])/num_lines #Total_G's in first position/Total_phasiRNA
            Norm_Count_phasiRNA[3][y] = float(seq_Comp_phasiRNA[3][y])/num_lines #Total_T's in first position/Total_phasiRNA  
 
#for non-phasiRNAs do the below
   
##################################################################################
##
##   Initialization of all variables
##   This routine counts freuncies of  A/C/G/T at each position
##
##################################################################################   

    seq_Comp_non_phasiRNA=np.zeros((4,n),int) #a 4x21 array of int zeros
    
    
    #read files
    input1=open(non_phasiRNA_sequences, 'r')
     
    for lines in input1:
        #print "Reading Line: %s" % (lines)
        sequence = lines.strip()
        if not lines: break

        for y in range (0,n):
            if (sequence [y] == 'A'):
                seq_Comp_non_phasiRNA [0][y] =  seq_Comp_non_phasiRNA [0][y]+1
            if (sequence [y] == 'C'):
                seq_Comp_non_phasiRNA [1][y] = seq_Comp_non_phasiRNA [1][y]+1
            if (sequence [y] == 'G'):
                seq_Comp_non_phasiRNA [2][y] = seq_Comp_non_phasiRNA [2][y]+1
            if (sequence [y] == 'T'):
                seq_Comp_non_phasiRNA [3][y] = seq_Comp_non_phasiRNA [3][y]+1                
 
########################################################################################
##
##   Initialization of all variables
##   This routine gets normalized frequency of nucleotides at each position in non-phasiRNAs.
##
########################################################################################    
           
    (rows,columns)=seq_Comp_non_phasiRNA.shape 
     
    Norm_Count_non_phasiRNA=np.zeros((4,columns),float) #a 4x21 array of int zeros
    
    for x in range(rows):
        for y in range(columns):
          
            Norm_Count_non_phasiRNA[0][y] = float(seq_Comp_non_phasiRNA[0][y])/num_lines1 #Total_A's in first position/Total_phasiRNA
            Norm_Count_non_phasiRNA[1][y] = float(seq_Comp_non_phasiRNA[1][y])/num_lines1 #Total_C's in first position/Total_phasiRNA
            Norm_Count_non_phasiRNA[2][y] = float(seq_Comp_non_phasiRNA[2][y])/num_lines1 #Total_G's in first position/Total_phasiRNA
            Norm_Count_non_phasiRNA[3][y] = float(seq_Comp_non_phasiRNA[3][y])/num_lines1 #Total_T's in first position/Total_phasiRNA       


####################################################################################################
##
##   This routine gets  frequency of nucleotides at each position in phasiRNAs and in non-phasiRNAs.
##   Creating 927x(4x21) or 927x(4x24) vector
##
####################################################################################################      
    
    # vector storing the information of frequencies of a1,c1,g1,t1,a2,c2,g2,t2..a21,c21,g21,t21 for phasiRNAs and non-phasiRNAs.    
    master_positional_vector_phasiRNA= np.zeros((num_lines,4*n),int)
    master_positional_vector_non_phasiRNA= np.zeros((num_lines1,4*n),int) 

    i=0
    input=open(phasiRNA_sequences, 'r')
    for lines in input:

        sequence = lines.strip()
        if not lines: break
        for y in range (0,n):
            if (sequence [y] == 'A'):
                master_positional_vector_phasiRNA[i][4*y]=1
            elif (sequence [y] == 'C'):
                master_positional_vector_phasiRNA[i][4*y+1]=1
            elif (sequence [y] == 'G'):
                master_positional_vector_phasiRNA[i][4*y+2]=1
            elif (sequence [y] == 'T'):
                master_positional_vector_phasiRNA[i][4*y+3]=1
            else:
                continue         
        i+=1
        
     #for non-phasiRNA         
    i=0
    input1=open(non_phasiRNA_sequences, 'r')
    for lines in input1:
            
        sequence = lines.strip()
        if not lines: break
        for y in range (0,n):
            if (sequence [y] == 'A'):      
                master_positional_vector_non_phasiRNA[i][4*y]=1
            elif (sequence [y] == 'C'):
                master_positional_vector_non_phasiRNA[i][4*y+1]=1
            elif (sequence [y] == 'G'):
                master_positional_vector_non_phasiRNA[i][4*y+2]=1
            elif (sequence [y] == 'T'):
                master_positional_vector_non_phasiRNA[i][4*y+3]=1
            else:
                continue          
        i+=1
        
        
        
        
    row,columns=master_positional_vector_phasiRNA.shape #926x84 vector
    p_value=[]
    
   
    alphabets=['A','C','G','T']
    for n in range (columns):
        group1=list(master_positional_vector_phasiRNA[:,n])
        group2=list(master_positional_vector_non_phasiRNA[:,n])
        p_val=RankSumTest(group1,group2)
        p_value.append(p_val)           
            
    
    i=1        
    significant_positions=[]
    print "significant_position_specific_base_usages at %s significance level\n" % (significance)
    for n in range(0,len(p_value)):
        
        
        
        if (p_value[n] < significance):        
            if(n%4==0):
                significant_positions.append(alphabets[0]+str(i))
                print alphabets[0],i,"-> p-val = ",p_value[n]
            elif(n%4==1):
                significant_positions.append(alphabets[1]+str(i))
                print alphabets[1],i,"-> p-val = ",p_value[n]
            elif(n%4==2):
                significant_positions.append(alphabets[2]+str(i))
                print alphabets[2],i,"-> p-val = ",p_value[n]
            elif(n%4==3):
                significant_positions.append(alphabets[3]+str(i))
                print alphabets[3],i,"-> p-val = ",p_value[n]
            else:
                continue
        if (n%4==3):
            i+=1
            
            

    significant_pos_freq_phasiRNA=[]
    significant_pos_freq_non_phasiRNA=[]
        
    for i in range(0,len(significant_positions)):
        
        pos=significant_positions[i] 
        if (pos[0]=="A"):
            pos_index=int(pos[1:])
            significant_pos_freq_phasiRNA.append(Norm_Count_phasiRNA[0,pos_index-1])
            significant_pos_freq_non_phasiRNA.append(Norm_Count_non_phasiRNA[0,pos_index-1])
        elif (pos[0]=="C"):
            pos_index=int(pos[1:])
            significant_pos_freq_phasiRNA.append(Norm_Count_phasiRNA[1,pos_index-1])
            significant_pos_freq_non_phasiRNA.append(Norm_Count_non_phasiRNA[1,pos_index-1])
        
        elif (pos[0]=="G"):
            pos_index=int(pos[1:])
            significant_pos_freq_phasiRNA.append(Norm_Count_phasiRNA[2,pos_index-1])
            significant_pos_freq_non_phasiRNA.append(Norm_Count_non_phasiRNA[2,pos_index-1])
            
        elif (pos[0]=="T"):
            pos_index=int(pos[1:])
            significant_pos_freq_phasiRNA.append(Norm_Count_phasiRNA[3,pos_index-1])
            significant_pos_freq_non_phasiRNA.append(Norm_Count_non_phasiRNA[3,pos_index-1])
        else:
            continue
            

    drawPFM(Norm_Count_phasiRNA,Norm_Count_non_phasiRNA,n,sorted(significant_positions))     
   
    input.close()
    input1.close()    
    


if __name__ == '__main__':   #calls the main function

    phasiRNA_sequences=positiveSet
    non_phasiRNA_sequences=negativeSet
    k_mer_list=[]
    n=size
    significance=p_val

    main(phasiRNA_sequences,non_phasiRNA_sequences,k_mer_list,n,significance)

    print "----------------------\nDone"
    sys.exit()
