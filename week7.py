'''Progressive failure analysis of Laminated Plates 
Two different supported panels are to be analyses 1) square panel a = b 2) rectangual panel. b = 2a
Load is uniform pressure p0 applied to the entire panel or a quarter of the panel. 
The deflection at the centre of the panel is calculated using classical plate theory. The moments at the centre and corners of the panel are also calculated.'''

import numpy as np
import functions as func




###INPUTS##############################################################################
# Input the material properties for the composite material
E1 = 181e3  # MPa
E2 = 10.3e3  # MPa
nu12 = 0.28
G12 = 7.17e3 # MPa

#change geometry of the plate
a = 1000 #mm
b = 1000 #mm
#a = 500 #mm

#chnage laminate layup and thickness
laminate = [0, 0,0,0,0,0,0,0]  # laminate layup in degrees
laminate = np.full_like(laminate, 90) # 8 plies of 0 degrees
laminate = [0,90,0,90,90,0,90,0] # 8 plies of 0 and 90 degrees
laminate = [45,-45,45,-45,-45,45,-45,45] # 8 plies of 45 and -45 degrees
#uniform thickness
t = 0.5 # 5thickness of each ply mm


#plate is fully loaded with uniform pressure p0 in MPa
p0 = 0.001 #MPa

alpha1 = 0
alpha2 = 1
beta1 = 0
beta2 = 1




#calculate D
Q = func.ply(E1, E2, nu12, G12) #calculate the local stiffness matrix for a single ply
A, B, D, z = func.laminate_properties(laminate, t, Q)

#calculate deflection at centre of plate (a/2, b/2) using orthropic plate theory

w,Mx,My,Mxy = func.orthropic_plate_arb(a, b, p0, D,11,11,alpha1,beta1,alpha2,beta2)

print(f"Deflection at centre of plate: {w} mm")
print(f"Moment Mx at centre of plate: {Mx} N")
print(f"Moment My at centre of plate: {My} N")
print(f"Moment Mxy at corner of plate: {Mxy} N ")


#only a quarter of the plate is loaded, load is same

alpha1 = 0.25
beta1 = 0.25
alpha2 = 0.75
beta2 = 0.75

w,Mx,My,Mxy = func.orthropic_plate_arb(a, b, p0, D,11,11,alpha1,beta1,alpha2,beta2)

print(f"Deflection at centre of plate: {w} mm")
print(f"Moment Mx at centre of plate: {Mx} N")
print(f"Moment My at centre of plate: {My} N")
print(f"Moment Mxy at corner of plate: {Mxy} N ")









 
    
   