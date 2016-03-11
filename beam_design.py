# -*- coding: utf-8 -*-
"""
Simple Span Beam Design

Author: Margaret Wang

This program is developed to design a beam for the geometry and loading provided by a user.
Design uses LRFD factors in calculating capacity and loading requirements.
The shape database used is limited to W shapes.

For this design, gravity loading in the major axis is considered.
Moment, shear, and axial load is considered.
+ gravity loading is downward
+ axial loading is compressive
"""

### Import libraries

# importing division from the 'future' release of Python (i.e. Python 3)
from __future__ import division

# importing libraries
import os
import pandas as pd
import math

### Define path to database of shapes for design

# defining path of directory where shapes database lives
# os.getcwd() # provides current working directory
path = "/Users/Margaret/Desktop/structural_data/"

os.chdir(os.path.dirname(path))
os.path.dirname(os.path.realpath('__file__'))

# importing database
shapes = pd.read_csv("AISC Shapes Database v14.1.csv")


### Define constants

# steel material properties
Fy = 50 #ksi
Fu = 65 #ksi
E = 29000 #ksi

# LRFD Factors
phi_yield = 0.90
phi_rupture = 0.75
phi_compression = 0.90
phi_flexure = 0.90
phi_shear = 1.0

#####
### Define Functions
#####

def check_value(num):
    ''' Checks if number entered is a float or not.
    '''
    valid = False    
    
    while not valid:
        try:
            float(num)
            valid = True
        except ValueError:
            num = raw_input("Sorry, not a number. Please input a number: ")
    
    value = float(num)
    return value

def check_string(char):
    ''' Checks if number entered is a float or not.
    '''
    valid = False
    
    while not valid:
        valid = char.upper() == 'Y' or char.upper() == 'N'
    
        if not valid:
            char = raw_input("Sorry, not a valid input. Please re-enter: ")
    
    return char.upper()

def request_areaload(load_case):
    ''' Requests area load input from user.
    '''
    area_load = raw_input("What is the " + load_case + " load in psf? ")
    area_value = check_value(area_load)
                    
    return area_value


def request_axialload(load_case):
    ''' Requests point load input from user.
    '''
    
    point_load = raw_input("What is the " + load_case + " axial load in kips? ")    
    point_value = check_value(point_load)
    
    return point_value


def request_geometry():
    ''' Requests geometry from user.
    '''
    
    length = raw_input("Enter length of beam in ft: ")
    length_value = check_value(length)
    
    spacing = raw_input("Enter beam spacing in ft: ")
    spacing_value = check_value(spacing)    
    
    unbraced_L = raw_input("Enter unbraced length of beam in ft (0 for fully braced): ")
    unbraced_L_value = check_value(unbraced_L)    
    
    K = raw_input("Enter effective length factor, K: ")
    K_value = check_value(K)    
    
    return length_value, spacing_value, unbraced_L_value, K_value   


def request_deflcriteria():
    ''' Requests total and live load deflection criteria from user.
    '''
    default = raw_input("Use default maximum total deflection of L/240 and "
                        " maximum live load deflection of L/360? \n Enter Y/N: ")
    default_YN = check_string(default)
    
    if default_YN == 'Y':
        total_defl_value = 240.0
        live_defl_value = 360.0
    else: 
        total_defl = raw_input("Enter maximum total deflection of beam, L/: ")
        total_defl_value = check_value(total_defl)    
        
        live_defl = raw_input("Enter maximum live load deflection of beam, L/: ")    
        live_defl_value = check_value(live_defl)
    
    return total_defl_value, live_defl_value

    
def input_loading():
    ''' Inputs gravity loading requirements.
    '''
    
    dead_area = request_areaload("dead")
    live_area = request_areaload("live")
    snow_area = request_areaload("snow")
    
    dead_axial = request_axialload("dead")
    live_axial = request_axialload("live")
    snow_axial = request_axialload("snow")
    
    return dead_area, live_area, snow_area, dead_axial, live_axial, snow_axial    


def max_loading(D,L,S):
    ''' ASCE 7-05 load combinations, focusing on gravity
    '''
    
    lc1 = 1.4*D
    lc2 = 1.2*D + 1.6*L + 0.5*S
    lc3 = 1.2*D + 1.0*L + 1.6*S
    
    strength = max(lc1, lc2, lc3)
    serviceability = D + L + S
    
    return strength, serviceability

def max_loading_df(df,col,D,L,S,spacing):
    ''' ASCE 7-05 load combinations, focusing on gravity
    '''
    
    df['lc1'] = 1.4*(D + df[col]/spacing)
    df['lc2'] = 1.2*(D + df[col]/spacing) + 1.6*L + 0.5*S
    df['lc3'] = 1.2*(D + df[col]/spacing) + 1.0*L + 1.6*S
    
    df['strength_load'] = df[['lc1', 'lc2', 'lc3']].max(axis=1)
    df['service_load'] = (D + df[col]/spacing) + L + S
    
    return df

   
def max_shaperatios(df, axial=False):
    ''' Checks compactness and slenderness of flanges and webs in
    flexure and in axial compression
    '''
    compact_class_flanges_flex = 0.38*math.sqrt(E/Fy)
    #noncompact_class_flanges_flex = 1.0*math.sqrt(E/Fy)
    df = df[df.bf_2tf < compact_class_flanges_flex]
    
    compact_class_webs_flex = 3.76*math.sqrt(E/Fy)
    df = df[df.d_2k_tw < compact_class_webs_flex]
    
    if axial:
        compact_class_flanges_comp = 0.56*math.sqrt(E/Fy)
        df = df[df.bf_2tf < compact_class_flanges_comp]
    
        compact_class_webs_comp = 1.49*math.sqrt(E/Fy)
        df = df[df.d_2k_tw < compact_class_webs_comp]  
    
    return df


def compare_cols(df,col1,col2):
    subset = df[df[col1] > df[col2]]
    return subset


####
## Design Process
####

### Request information from users

D_psf, L_psf, S_psf, D_kip, L_kip, S_kip = input_loading() # generates loads for load cases
L, spacing, Lb, K = request_geometry() # generates spacing


### Prepare database of shapes to check
# isolating W shapes
W_shapes = shapes[shapes['Type']=='W']

# limiting to necessary elements
W_shapes = W_shapes[['AISC_Manual_Label','W','A','d','bf','tw','tf','kdes','bf/2tf','h/tw',
                    'Ix','Zx','Sx','rx','Iy','Zy','Sy','ry','J','Cw','rts','ho']]
W_shapes.rename(columns = {'AISC_Manual_Label':'Shape','bf/2tf':'bf_2tf','h/tw':'h_tw'}, inplace=True) # renames columns
#W_shapes = W_shapes.convert_objects(convert_numeric=True) # convert numeric (deprecated)
W_shapes = W_shapes.apply(pd.to_numeric, errors='ignore')  # convert numeric


### Perform design

### develop design requirements
strength_psf, service_psf = max_loading(D_psf,L_psf,S_psf) # generates area load before beam laod
strength_kip, service_kip = max_loading(D_kip, L_kip, S_kip) # generates axial point loading

strength_w = strength_psf*spacing/1000
service_w = service_psf*spacing/1000
live_w = L_psf*spacing/1000

W_shapes = max_loading_df(W_shapes,'W',D_psf,L_psf,S_psf,spacing) # generates line loading
W_shapes['strength_w'] = W_shapes['strength_load']*spacing/1000
W_shapes['service_w'] = W_shapes['service_load']*spacing/1000

strength_kip, service_kip = max_loading(D_kip, L_kip, S_kip) #generates axial point loading
W_shapes['strength_k'] = strength_kip
W_shapes['service_k'] = service_kip

W_shapes['moment_ult'] = W_shapes.strength_w*L**2/8
W_shapes['shear_ult'] = W_shapes.strength_w*L/2
W_shapes['axial_ult'] = strength_kip

W_shapes['moment_nom'] = W_shapes.service_w*L**2/8
W_shapes['shear_nom'] = W_shapes.service_w*L/2
W_shapes['axial_nom'] = service_kip


moment_ult = strength_w*L**2/8
shear_ult = strength_w*L/2
axial_ult = strength_kip

moment_nom = service_w*L**2/8
shear_nom = service_w*L/2
axial_nom = service_kip


print('The LRFD design values (without beam weight) are: \n' +  
'Moment:  %0.2f kip-ft\n' +
'Shear: %0.2f kips\n' +
'Axial Load: %0.2f kips\n') % (moment_ult, shear_ult, axial_ult)

print('The unfactored design values (without beam weight) are: \n' +  
'Moment:  %0.2f kip-ft\n' +
'Shear: %0.2f kips\n' +
'Axial Load: %0.2f kips\n')  % (moment_nom, shear_nom, axial_nom)

W_shapes['d_2k_tw'] = (W_shapes.d - 2*W_shapes.kdes) / W_shapes.tw # create (d-2*k)/tw column
W_shapes['A_web'] = W_shapes.d*W_shapes.tw # create area of web column = d*tw

# defines subset of shapes that fit shape ratio requirements in flexure and compression
axial = service_kip > 0
W_shapes = max_shaperatios(W_shapes, axial) 

####
## Deflection Design
####

### ∆ = 5wL^4/384EI
### Imin = 5wL^4/384EI∆

total_defl_crit, live_defl_crit = request_deflcriteria()

total_defl_limit = L*12/total_defl_crit
live_defl_limit = L*12/live_defl_crit

print('The deflection criteria translates to: \n' +  
'Total (L/%d):  %0.2f in\n' +
'Live (L/%d): %0.2f in\n')  % (total_defl_crit, total_defl_limit, live_defl_crit, live_defl_limit)


# Calculate minimum Ix
W_shapes['Imin'] = 5*(W_shapes.service_w/12)*(L*12)**4/(384*E*(L*12/total_defl_crit))
Imin_L = 5*(live_w/12)*(L*12)**4/(384*E*(L*12/live_defl_crit))

# remove from list sections where the live load deflection governs
W_shapes.loc[W_shapes.Imin < Imin_L, 'Imin'] = Imin_L

# compares shape Ix to Imin
W_shapes = compare_cols(W_shapes,'Ix','Imin')

# Deflections
W_shapes['total_defl'] = 5*(W_shapes.service_w/12)*(L*12)**4/(384*E*W_shapes.Ix)
W_shapes['live_defl'] = 5*(live_w/12)*(L*12)**4/(384*E*W_shapes.Ix)

# Return database of shapes that satisfy criteria
valid_defl_shapes = W_shapes[W_shapes.Ix > W_shapes.Imin].sort_values(by=['Ix']) # all valid shapes
top_10_defl_shapes = valid_defl_shapes[['Shape','Ix', 'total_defl', 'live_defl']].head(10)
lightest_defl_shape = W_shapes.loc[W_shapes[W_shapes.Ix > W_shapes.Imin].W.idxmin(),'Shape'] # smallest shape

print 'Lightest shape for deflection is: ' + lightest_defl_shape
print 'Top 10 efficient shapes for deflection: \n' + top_10_defl_shapes.to_string(index=False)

####
## Flexural Design (Chapter C and F)
####

### Amplified first order elastic analysis
Cm = 1.0 # coefficient assuming no lateral translation, Cm = 0.6 - 0.4 M1/M2
alpha = 1.0 # for LRFD

Pe1 = math.pi**2*E*W_shapes.Ix/(K*L*12)**2
W_shapes['B1'] = Cm/(1-alpha*W_shapes.axial_ult/Pe1)
W_shapes.loc[W_shapes['B1']<1, 'B1'] = 1

W_shapes['moment_ult'] = W_shapes.B1*W_shapes.moment_ult

### Flexural yielding
### Mn = FyZx
### phiMn = phi_flexure*Fy*Zx
W_shapes['M_fy'] = phi_flexure*Fy*W_shapes.Zx/12

### Lateral-torsional buckling
c = 1.0
Lb_in = Lb*12

Cb = 1.0
h0 = W_shapes.d-W_shapes.tf

W_shapes['Lp'] = 1.76*W_shapes.ry*(E/Fy)**0.5
W_shapes['Lr'] = 1.95*W_shapes.rts*E/(0.7*Fy)*(W_shapes.J*c/(W_shapes.Sx*h0)+((W_shapes.J*c/(W_shapes.Sx*h0))**2+6.76*(0.7*Fy/E)**2)**0.5)**0.5

## Case: Lb < Lp, no lateral torsional buckling, yield moment governs
no_LTB = W_shapes['Lp'] > Lb_in # boolean mask for Lp > Lb
#W_shapes.loc[no_LTB, 'Shape'] # list of shapes that encounter lateral torsional buckling
W_shapes.loc[no_LTB, 'M_LTB'] = W_shapes['M_fy'] #dataframe of no_LTB shapes

## Case: Lp < Lb < Lr, inelastic lateral torsional buckling
LTB_ie = (W_shapes['Lp'] < Lb_in) & (W_shapes['Lr'] > Lb_in) # boolean mask
M_LTB_ie = Cb*(W_shapes.M_fy - (W_shapes.M_fy-0.7*Fy*W_shapes.Sx)*(Lb_in - W_shapes.Lp)/(W_shapes.Lr-W_shapes.Lp))
W_shapes.loc[LTB_ie, 'M_LTB'] = phi_flexure*M_LTB_ie/12 # assign to frame

## Case: Lb > Lr, elastic lateral torsional buckling
LTB_e = W_shapes['Lr'] < Lb_in # boolean mask
Fcr = Cb*math.pi**2*E/(Lb_in/W_shapes.rts)**2*(1+0.78*W_shapes.J*c/(W_shapes.Sx*h0)*(Lb_in/W_shapes.rts)**2)**0.5
M_LTB_e = Fcr*W_shapes.Sx
W_shapes.loc[LTB_e, 'M_LTB'] = phi_flexure*M_LTB_e/12

W_shapes['M'] = W_shapes[['M_fy', 'M_LTB']].min(axis=1)

# available shapes that fit design criteria

valid_flex_shapes = W_shapes[W_shapes.M > W_shapes.moment_ult].sort_values(by=['Zx']) # all valid shapes
top_flex_shapes = valid_flex_shapes[['Shape','Zx','M']].head(10)
lightest_flex_shape = W_shapes.loc[W_shapes[W_shapes.M > W_shapes.moment_ult].W.idxmin(),'Shape'] # smallest shape

print 'Lightest shape for flexure is: ' + lightest_flex_shape
print 'Top 10 efficient shapes for flexure: \n' + top_flex_shapes.to_string(index=False)

####
## Shear Design (Chapter G)
####

## Vn = 0.6*Fy*Aw*Cv
## phiVn = phi_shear*0.6*Fy*Aw*Cv
Cv = 1.0 # From user note Chapter G2.1
W_shapes['V'] = phi_shear*0.6*Fy*Cv*W_shapes.A_web

# available shapes that fit design criteria
valid_shear_shapes = W_shapes[W_shapes.V > W_shapes.shear_ult].sort_values(by=['A_web']) # all valid shapes
top_shear_shapes = valid_shear_shapes[['Shape','A_web','V']].head(10)
lightest_shear_shape = W_shapes.loc[W_shapes[W_shapes.V > W_shapes.shear_ult].W.idxmin(),'Shape'] # smallest shape

print 'Lightest shape for shear is: ' + lightest_shear_shape
print 'Top 10 efficient shapes for shear: \n' + top_shear_shapes.to_string(index=False)


####
## Axial Loading in Compression (Chapter E)
####

# Minimum Slenderness
KL_r200 = K*L*12/W_shapes.rx < 200

# remove from list sections where the live load deflection governs
W_shapes = W_shapes.loc[KL_r200]

# Major Axis Buckling
W_shapes['Fe'] = math.pi**2*E/(K*L*12/W_shapes.rx)**2 # elastic buckling stress
KL_r = 4.71*(E/Fy)**0.5 # slenderness
W_shapes['Fc'] = 0.877*W_shapes.Fe # critical stress

## Case: KL/R < KL_r
KL_r_crit = K*L*12/W_shapes.rx < KL_r # boolean
Fc_crit = Fy*0.658**(Fy/W_shapes.Fe)
W_shapes.loc[KL_r_crit, 'Fc'] = Fc_crit # replacing other critical stress

W_shapes['P'] = phi_compression*W_shapes.Fc*W_shapes.A

####
## Combined Forces (Chapter H)
####

W_shapes['P_M'] = W_shapes.axial_ult/(2*W_shapes.P) + W_shapes.moment_ult/W_shapes.M
combined_forces = W_shapes.P_M < 1

W_shapes[combined_forces]

valid_combined_shapes = W_shapes[combined_forces].sort_values(by=['A']) # all valid shapes
top_combined_shapes = valid_combined_shapes[['Shape','Ix','Zx','A','M','V','P']].head(10)

lightest_combined_shape = W_shapes.loc[W_shapes[combined_forces].W.idxmin(),'Shape'] # smallest shape

print 'Lightest shape for combined forces is: ' + lightest_combined_shape
print 'Top 10 lightest shapes for combined forces: \n' + top_combined_shapes.to_string(index=False)

print 'The selected shape for the loading: '
print W_shapes[W_shapes.Shape.isin([lightest_combined_shape])][['Shape', 
               'Ix','total_defl','live_defl','Zx','moment_ult','M','shear_ult','V','axial_ult','P', 'P_M']]
              