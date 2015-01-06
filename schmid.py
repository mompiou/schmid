# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 17:05:02 2014

@author: mompiou
"""
from __future__ import division
import numpy as np
from Tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib import pyplot as plt
import ttk
import sys
from tkFileDialog import *
import os
#import pdb


pi=np.pi


def d(pole1,pole2,pole3):
	global M,V,D,Dstar
	
	a=eval(a_entry.get())
	b=eval(b_entry.get())
	c=eval(c_entry.get())
	alp=eval(alp_entry.get())
	bet=eval(bet_entry.get())
	gam=eval(gam_entry.get())
	alp=alp*pi/180;
	bet=bet*pi/180;
	gam=gam*pi/180;
	V=a*b*c*np.sqrt(1-(np.cos(alp)**2)-(np.cos(bet))**2-(np.cos(gam))**2+2*b*c*np.cos(alp)*np.cos(bet)*np.cos(gam))
	D=np.array([[a,b*np.cos(gam),c*np.cos(bet)],[0,b*np.sin(gam),c*(np.cos(alp)-np.cos(bet)*np.cos(gam))/np.sin(gam)],[0,0,V/(a*b*np.sin(gam))]])
	Dstar=np.transpose(np.linalg.inv(D))
	G=np.array([[a**2,a*b*np.cos(gam),a*c*np.cos(bet)],[a*b*np.cos(gam),b**2,b*c*np.cos(alp)],[a*c*np.cos(bet),b*c*np.cos(alp),c**2]])    
	ds=(np.sqrt(np.dot(np.array([pole1,pole2,pole3]),np.dot(np.linalg.inv(G),np.array([pole1,pole2,pole3])))))
	return ds



    
def pole(pole1,pole2,pole3):
	global M,V,D,Dstar,G, alp, bet, gam    
	
	a=eval(a_entry.get())
	b=eval(b_entry.get())
	c=eval(c_entry.get())
	alp=eval(alp_entry.get())
	bet=eval(bet_entry.get())
	gam=eval(gam_entry.get())
	alp=alp*pi/180;
	bet=bet*pi/180;
	gam=gam*pi/180;
	G=np.array([[a**2,a*b*np.cos(gam),a*c*np.cos(bet)],[a*b*np.cos(gam),b**2,b*c*np.cos(alp)],[a*c*np.cos(bet),b*c*np.cos(alp),c**2]])     
	v=d(pole1,pole2,pole3)
	N=np.array([pole1,pole2,pole3])
	if np.abs(alp-pi/2)<0.001 and np.abs(bet-pi/2)<0.001 and np.abs(gam-2*pi/3)<0.001:
		N=np.array([[pole1,pole2,pole3],[pole1,pole2,-pole3],[pole2,pole1,pole3],[pole2,pole1,-pole3],[-pole1-pole2,pole2,pole3],[-pole1-pole2,pole2,-pole3],[pole1,-pole1-pole2,pole3],[pole1,-pole1-pole2,-pole3],[pole2,-pole1-pole2,pole3],[pole2,-pole1-pole2,-pole3],[-pole1-pole2,pole1,pole3],[-pole1-pole2,pole1,-pole3]])
	else:	
        	if np.abs(d(pole1,pole2,-pole3)-v)<0.001:
                    N=np.vstack((N,np.array([pole1,pole2,-pole3])))
        	if np.abs(d(pole1,-pole2,pole3)-v)<0.001:
        			N=np.vstack((N,np.array([pole1,-pole2,pole3])))
        	if np.abs(d(-pole1,pole2,pole3)-v)<0.001:
        			N=np.vstack((N,np.array([-pole1,pole2,pole3])))
        	if np.abs(d(pole2,pole1,pole3)-v)<0.001:
        			N=np.vstack((N,np.array([pole2,pole1,pole3])))
        	if np.abs(d(pole2,pole1,-pole3)-v)<0.001:
        		N=np.vstack((N,np.array([pole2,pole1,-pole3])))
        	if np.abs(d(pole2,-pole1,pole3)-v)<0.001:
        		N=np.vstack((N,np.array([pole2,-pole1,pole3])))
        	if np.abs(d(-pole2,pole1,pole3)-v)<0.001:
        		N=np.vstack((N,np.array([-pole2,pole1,pole3])))
        	if np.abs(d(pole2,pole3,pole1)-v)<0.001:
        		N=np.vstack((N,np.array([pole2,pole3,pole1])))
        	if np.abs(d(pole2,pole3,pole1)-v)<0.001:
        		N=np.vstack((N,np.array([pole2,pole3,-pole1])))
        	if np.abs(d(pole2,-pole3,pole1)-v)<0.001:
        		N=np.vstack((N,np.array([pole2,-pole3,pole1])))
        	if np.abs(d(-pole2,pole3,pole1)-v)<0.001:
        		N=np.vstack((N,np.array([-pole2,pole3,pole1])))
        	if np.abs(d(pole1,pole3,pole2)-v)<0.001:
        		N=np.vstack((N,np.array([pole1,pole3,pole2])))
        	if np.abs(d(pole1,pole3,-pole2)-v)<0.001:
        		N=np.vstack((N,np.array([pole1,pole3,-pole2])))
        	if np.abs(d(pole1,-pole3,pole2)-v)<0.001:
        		N=np.vstack((N,np.array([pole1,-pole3,pole2])))
        	if np.abs(d(-pole1,pole3,pole2)-v)<0.001:
        		N=np.vstack((N,np.array([-pole1,pole3,pole2])))
        	if np.abs(d(pole3,pole1,pole2)-v)<0.001:
        		N=np.vstack((N,np.array([pole3,pole1,pole2])))
        	if np.abs(d(pole3,pole1,-pole2)-v)<0.001:
        		N=np.vstack((N,np.array([pole3,pole1,-pole2])))
        	if np.abs(d(pole3,-pole1,pole2)-v)<0.001:
        		N=np.vstack((N,np.array([pole3,-pole1,pole2])))
        	if np.abs(d(-pole3,pole1,pole2)-v)<0.001:
        		N=np.vstack((N,np.array([-pole3,pole1,pole2])))
        	if np.abs(d(pole3,pole2,pole1)-v)<0.001:
        		N=np.vstack((N,np.array([pole3,pole2,pole1])))
        	if np.abs(d(pole3,pole2,-pole1)-v)<0.001:
        		N=np.vstack((N,np.array([pole3,pole2,-pole1])))
        	if np.abs(d(pole3,-pole2,pole1)-v)<0.001:
        		N=np.vstack((N,np.array([pole3,-pole2,pole1])))
        	if np.abs(d(pole3,pole2,pole1)-v)<0.001:
        		N=np.vstack((N,np.array([pole3,pole2,pole1])))
	
		  
	return N


def rotation(phi1,phi,phi2):
   phi1=phi1*pi/180;
   phi=phi*pi/180;
   phi2=phi2*pi/180;
   R=np.array([[np.cos(phi1)*np.cos(phi2)-np.cos(phi)*np.sin(phi1)*np.sin(phi2),
            -np.cos(phi)*np.cos(phi2)*np.sin(phi1)-np.cos(phi1)*
            np.sin(phi2),np.sin(phi)*np.sin(phi1)],[np.cos(phi2)*np.sin(phi1)
            +np.cos(phi)*np.cos(phi1)*np.sin(phi2),np.cos(phi)*np.cos(phi1)
            *np.cos(phi2)-np.sin(phi1)*np.sin(phi2), -np.cos(phi1)*np.sin(phi)],
            [np.sin(phi)*np.sin(phi2), np.cos(phi2)*np.sin(phi), np.cos(phi)]],float)
   return R


def schmid(b,n):
    global D, Dstar,M, alp, bet, gam  
    if np.abs(alp-pi/2)<0.001 and np.abs(bet-pi/2)<0.001 and np.abs(gam-2*pi/3)<0.001:
        b2=np.array([0,0,0])        
        b2[0]=2*b[0]+b[1]
        b2[1]=2*b[1]+b[0]
        b2[2]=b[2]
        bpr=np.dot(D, b2)
    else:
        bpr=np.dot(D, b)
     
    npr=np.dot(Dstar,n)
   
    npr2=np.dot(M,npr)
    bpr2=np.dot(M,bpr)
    t0=eval(t0_entry.get())
    t1=eval(t1_entry.get())
    t2=eval(t2_entry.get())
    T=np.array([t0,t1,t2])/np.linalg.norm(np.array([t0,t1,t2]))
    anglen=np.arccos(np.dot(npr2,T)/np.linalg.norm(npr2))
    angleb=np.arccos(np.dot(bpr2,T)/np.linalg.norm(bpr2))
    
    s=np.cos(anglen)*np.cos(angleb)

    return s

def prod_scal(c1,c2):
    global M, Dstar, D, alp, bet,gam
    if np.abs(alp-pi/2)<0.001 and np.abs(bet-pi/2)<0.001 and np.abs(gam-2*pi/3)<0.001:
        c2p=np.array([0,0,0])        
        c2p[0]=2*c2[0]+c2[1]
        c2p[1]=2*c2[1]+c2[0]
        c2p[2]=c2[2]
        c2c=np.dot(D, c2p)
    else:
        c2c=np.dot(D, c2)
    
    c1c=np.dot(Dstar,c1)
    
    
    p=np.dot(c1c,c2c)
    
    #print(p)
    return p

def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))
    

def main():
	global M
        u=eval(b1_entry.get())   
        v=eval(b2_entry.get())
        w=eval(b3_entry.get())
        h=eval(n1_entry.get())   
        k=eval(n2_entry.get())
        l=eval(n3_entry.get())
        phi1=eval(phi1_entry.get())   
        phi=eval(phi_entry.get())
        phi2=eval(phi2_entry.get())
        	
        	
        B=pole(u,v,w)
        N=pole(h,k,l)
        #print(B,N)
        P=np.array([0,0,0,0,0,0,0])
        M=rotation(phi1,phi,phi2)
        Listbox_schmid.delete(0,END) 
     
	for i in range(0,np.shape(N)[0]):
		for j in range(0,np.shape(B)[0]):
                  if np.abs(prod_scal(N[i,:],B[j,:]))<0.0001:
                     s=schmid(B[j,:],N[i,:])
#                     if s>0:
                     R=np.array([s,N[i,0],N[i,1],N[i,2],B[j,0],B[j,1],B[j,2]])
                     P=np.vstack((P,R))
                     
        P=np.delete(P, (0), axis=0)
        P=unique_rows(P)
        P=-P
        P.view('float64,i8,i8,i8,i8,i8,i8').sort(order=['f0'], axis=0)
        P=-P
        #print(P)
        Listbox_schmid.insert(END,( 's', '|' ,'n', '|', 'b' ))
        for k in range(0,np.shape(P)[0]):                                   
            Listbox_schmid.insert(END,(np.around(P[k,0],decimals=3), '|' ,np.int(P[k,1]), np.int(P[k,2]),np.int(P[k,3]), '|', np.int(P[k,4]), np.int(P[k,5]),np.int(P[k,6])  ))
        


root = Tk()
root.wm_title("Schmid factor")
root.geometry('600x600+10+40')
root.configure(bg = '#BDBDBD')
#root.resizable(0,0)
#s=ttk.Style()
#s.theme_use('clam')
style = ttk.Style()
theme = style.theme_use()
default = style.lookup(theme, 'background')


a_cristal_label = Label (master=root)
a_cristal_label.place(relx=0.1,rely=0.06,height=19,width=12)
a_cristal_label.configure(text='''a''')

b_cristal_label = Label (master=root)
b_cristal_label.place(relx=0.1,rely=0.12,height=19,width=12)
b_cristal_label.configure(activebackground="#f9f9f9")
b_cristal_label.configure(activeforeground="black")
b_cristal_label.configure(foreground="black")
b_cristal_label.configure(highlightcolor="black")
b_cristal_label.configure(text='''b''')

c_cristal_label = Label (master=root)
c_cristal_label.place(relx=0.1,rely=0.18,height=19,width=11)
c_cristal_label.configure(activebackground="#f9f9f9")
c_cristal_label.configure(activeforeground="black")
c_cristal_label.configure(foreground="black")
c_cristal_label.configure(highlightcolor="black")
c_cristal_label.configure(text='''c''')

alp_cristal_label = Label (master=root)
alp_cristal_label.place(relx=0.06,rely=0.24,height=19,width=45)
alp_cristal_label.configure(activebackground="#f9f9f9")
alp_cristal_label.configure(activeforeground="black")
alp_cristal_label.configure(foreground="black")
alp_cristal_label.configure(highlightcolor="black")
alp_cristal_label.configure(text='''alpha''')

bet_cristal_label = Label (master=root)
bet_cristal_label.place(relx=0.06,rely=0.30,height=19,width=45)
bet_cristal_label.configure(activebackground="#f9f9f9")
bet_cristal_label.configure(activeforeground="black")
bet_cristal_label.configure(foreground="black")
bet_cristal_label.configure(highlightcolor="black")
bet_cristal_label.configure(text='''beta''')

gam_cristal_label = Label (master=root)
gam_cristal_label.place(relx=0.06,rely=0.36,height=19,width=50)
gam_cristal_label.configure(activebackground="#f9f9f9")
gam_cristal_label.configure(activeforeground="black")
gam_cristal_label.configure(foreground="black")
gam_cristal_label.configure(highlightcolor="black")
gam_cristal_label.configure(text='''gamma''')

a_entry = Entry (master=root)
a_entry.place(relx=0.15,rely=0.06,relheight=0.03,relwidth=0.06)
a_entry.configure(background="white")
a_entry.configure(insertbackground="black")

b_entry = Entry (master=root)
b_entry.place(relx=0.15,rely=0.12,relheight=0.03,relwidth=0.06)
b_entry.configure(background="white")
b_entry.configure(foreground="black")
b_entry.configure(highlightcolor="black")
b_entry.configure(insertbackground="black")
b_entry.configure(selectbackground="#c4c4c4")
b_entry.configure(selectforeground="black")

c_entry = Entry (master=root)
c_entry.place(relx=0.15,rely=0.18,relheight=0.03,relwidth=0.06)
c_entry.configure(background="white")
c_entry.configure(foreground="black")
c_entry.configure(highlightcolor="black")
c_entry.configure(insertbackground="black")
c_entry.configure(selectbackground="#c4c4c4")
c_entry.configure(selectforeground="black")

alp_entry = Entry (master=root)
alp_entry.place(relx=0.15,rely=0.24,relheight=0.03,relwidth=0.06)
alp_entry.configure(background="white")
alp_entry.configure(foreground="black")
alp_entry.configure(highlightcolor="black")
alp_entry.configure(insertbackground="black")
alp_entry.configure(selectbackground="#c4c4c4")
alp_entry.configure(selectforeground="black")

bet_entry = Entry (master=root)
bet_entry.place(relx=0.15,rely=0.30,relheight=0.03,relwidth=0.06)
bet_entry.configure(background="white")
bet_entry.configure(foreground="black")
bet_entry.configure(highlightcolor="black")
bet_entry.configure(insertbackground="black")
bet_entry.configure(selectbackground="#c4c4c4")
bet_entry.configure(selectforeground="black")

gam_entry = Entry (master=root)
gam_entry.place(relx=0.15,rely=0.36,relheight=0.03,relwidth=0.06)
gam_entry.configure(background="white")
gam_entry.configure(foreground="black")
gam_entry.configure(highlightcolor="black")
gam_entry.configure(insertbackground="black")
gam_entry.configure(selectbackground="#c4c4c4")
gam_entry.configure(selectforeground="black")

T_label = Label (master=root)
T_label.place(relx=0.08,rely=0.50,height=19,width=30)
T_label.configure(activebackground="#f9f9f9")
T_label.configure(activeforeground="black")
T_label.configure(foreground="black")
T_label.configure(highlightcolor="black")
T_label.configure(text='''T''')

t0_entry = Entry (master=root)
t0_entry.place(relx=0.15,rely=0.44,relheight=0.03,relwidth=0.06)
t0_entry.configure(background="white")
t0_entry.configure(foreground="black")
t0_entry.configure(highlightcolor="black")
t0_entry.configure(insertbackground="black")
t0_entry.configure(selectbackground="#c4c4c4")
t0_entry.configure(selectforeground="black")

t1_entry = Entry (master=root)
t1_entry.place(relx=0.15,rely=0.50,relheight=0.03,relwidth=0.06)
t1_entry.configure(background="white")
t1_entry.configure(foreground="black")
t1_entry.configure(highlightcolor="black")
t1_entry.configure(insertbackground="black")
t1_entry.configure(selectbackground="#c4c4c4")
t1_entry.configure(selectforeground="black")

t2_entry = Entry (master=root)
t2_entry.place(relx=0.15,rely=0.56,relheight=0.03,relwidth=0.06)
t2_entry.configure(background="white")
t2_entry.configure(foreground="black")
t2_entry.configure(highlightcolor="black")
t2_entry.configure(insertbackground="black")
t2_entry.configure(selectbackground="#c4c4c4")
t2_entry.configure(selectforeground="black")


b1_entry = Entry (master=root)
b1_entry.place(relx=0.25,rely=0.1,relheight=0.04,relwidth=0.08)
b1_entry.configure(background="white")

b2_entry = Entry (master=root)
b2_entry.place(relx=0.25,rely=0.15,relheight=0.04,relwidth=0.08)
b2_entry.configure(background="white")

b3_entry = Entry (master=root)
b3_entry.place(relx=0.25,rely=0.2,relheight=0.04,relwidth=0.08)
b3_entry.configure(background="white")

n1_entry = Entry (master=root)
n1_entry.place(relx=0.35,rely=0.1,relheight=0.04,relwidth=0.08)
n1_entry.configure(background="white")

n2_entry = Entry (master=root)
n2_entry.place(relx=0.35,rely=0.15,relheight=0.04,relwidth=0.08)
n2_entry.configure(background="white")

n3_entry = Entry (master=root)
n3_entry.place(relx=0.35,rely=0.2,relheight=0.04,relwidth=0.08)
n3_entry.configure(background="white")

phi1_entry = Entry (master=root)
phi1_entry.place(relx=0.5,rely=0.1,relheight=0.03,relwidth=0.07)
phi1_entry.configure(background="white")
phi1_entry.configure(foreground="black")
phi1_entry.configure(highlightcolor="black")
phi1_entry.configure(insertbackground="black")
phi1_entry.configure(selectbackground="#c4c4c4")
phi1_entry.configure(selectforeground="black")


phi_entry = Entry (master=root)
phi_entry.place(relx=0.5,rely=0.15,relheight=0.03,relwidth=0.07)
phi_entry.configure(background="white")
phi_entry.configure(foreground="black")
phi_entry.configure(highlightcolor="black")
phi_entry.configure(insertbackground="black")
phi_entry.configure(selectbackground="#c4c4c4")
phi_entry.configure(selectforeground="black")

label_euler = Label (master=root)
label_euler.place(relx=0.5,rely=0.05,height=19,width=70)
label_euler.configure(activebackground="#cccccc")
label_euler.configure(activeforeground="black")
label_euler.configure(foreground="black")
label_euler.configure(highlightcolor="black")
label_euler.configure(text='''Euler angle''')

phi2_entry = Entry (master=root)
phi2_entry.place(relx=0.5,rely=0.2,relheight=0.03,relwidth=0.07)
phi2_entry.configure(background="white")
phi2_entry.configure(foreground="black")
phi2_entry.configure(highlightcolor="black")
phi2_entry.configure(insertbackground="black")
phi2_entry.configure(selectbackground="#c4c4c4")
phi2_entry.configure(selectforeground="black")

schmid_button = Button (master=root)
schmid_button.place(relx=0.35,rely=0.3,relheight=0.04,relwidth=0.2)
schmid_button.configure(command=main)
schmid_button.configure(pady="0")
schmid_button.configure(text='''Schmid factor''')



Listbox_schmid = Listbox (master=root)
Listbox_schmid.place(relx=0.35,rely=0.35,height=190,width=256)
Listbox_schmid.configure(background="white")

burgers_label = Label (master=root)
burgers_label.place(relx=0.25,rely=0.05,height=19,width=60)
burgers_label.configure(text=''' b''')

normal_label = Label (master=root)
normal_label.place(relx=0.35,rely=0.05,height=19,width=60)
normal_label.configure(text=''' n''')


t0_entry.insert(1,0)
t1_entry.insert(1,1)
t2_entry.insert(1,0)
a_entry.insert(1,1)
b_entry.insert(1,1)
c_entry.insert(1,1)
alp_entry.insert(1,90)
bet_entry.insert(1,90)
gam_entry.insert(1,90)
b1_entry.insert(1,1)
b2_entry.insert(1,1)
b3_entry.insert(1,0)
n1_entry.insert(1,1)
n2_entry.insert(1,0)
phi1_entry.insert(1,0)
phi_entry.insert(1,0)
phi2_entry.insert(1,0)
n3_entry.insert(1,0)

mainloop()   

