import math as mh # R added
import numpy as np # R added

def single_harmonic(obs,time, h=1):
    "where obs is your data series"
    "Estimation of the amplitude and phase of a single harmonic"
    "Use in simple scenarios e.g sinusiodal annual cycle of temperature"
    "returns phase angle and amplitude"
    "See Wilks 2006 stats in atmos science for method"
    # R - equations for thsi are on p. 376
    pi=mh.pi
    cos=mh.cos
    sin=mh.sin
    atan=mh.atan
    mean=np.mean(obs)
    nobs=len(obs)
    list_A=[]
    list_B=[]
    for t,y in zip(time,obs):
        A=y*(cos((2*pi*t*h)/nobs))
        B=y*(sin((2*pi*t*h)/nobs))
        list_A.append(A)
        list_B.append(B)

    A_total=np.sum(list_A)
    B_total=np.sum(list_B)
    Ax=(2*(A_total))/nobs
    Bx=(2*(B_total))/nobs

    #defining C as the amplitude
    C=np.sqrt((Ax**2 + Bx**2))

    #phi as phase-angle
    print Ax
    if Ax>0:
        phi=atan((Bx/Ax))
    elif Ax<0:
        phi=atan((Bx/Ax)) #wasnt specified so I added this from above - I think it's correct
        if phi<pi:
            phi=phi+pi
        elif phi>pi:
            phi=phi-pi
    elif Ax==0:
        phi=pi/2

    phi=(phi/h)
    #convert phase angle in radians to degrees.
    ps_deg=mh.degrees(phi)

    return phi,ps_deg,C

def higher_harmonics(obs, time):

    "See description for simple_harmonic function. Here we include n/2 harmonics"
    "Amplitude variance explained by each harmonic is also computed"
    "returns list of amplitude variance explained, phase angles and amplitudes for each harmonic"
    pi=mh.pi
    cos=mh.cos
    sin=mh.sin
    atan=mh.atan
    mean=np.mean(obs)
    nobs=len(obs)
    harmonic_number=(nobs/2)+1
    S=np.var(obs)

    harm_funcs=np.arange(1,harmonic_number+1) # R edited, previously was (1,nobs)
    #this section works out the amplitude and phaseshift
    #1st work out A1 and B1 using formulae
    list_ps=[]
    list_C=[]
    for h in harm_funcs:
        list_A=[]
        list_B=[]
        for t,y in zip(time,obs):
            A=y*(cos((2*pi*h*t)/nobs))
            B=y*(sin((2*pi*h*t)/nobs))
            list_A.append(A)
            list_B.append(B)
        A_total=np.sum(list_A)
        B_total=np.sum(list_B)
        Ax=(2*(A_total))/nobs
        Bx=(2*(B_total))/nobs

        #defining C as the amplitude
        C=np.sqrt((Ax**2 + Bx**2))

        if Ax>0:
            phi=atan(Bx/Ax)
        elif Ax<0:
            if 0 < atan(Bx/Ax)+pi < 2*pi:
                phi=atan(Bx/Ax)+pi
            else:
                phi=atan(Bx/Ax)-pi
        elif Ax==0:
            phi=pi/2

        #convert phase angle in radians to degrees.
        ps_deg=mh.degrees(phi)
        list_ps.append(ps_deg)
        list_C.append(C)
    sum_amps=np.sum(list_C)
    ex_var_list=[]
    amp_var_list=[]


    for h in harm_funcs:
        index=h-1
        x=list_C[index]
        amp_var=(x/sum_amps)*100
        amp_var_list.append(amp_var)
        ex_variance=(((nobs/2)*(x**2))/((nobs-1)*S))*100

        ex_var_list.append(ex_variance)

    return list_C, list_ps, ex_var_list, amp_var_list


def higher_harmonics_fx(obs, time,nh=20):

    "See description for simple_harmonic function"
    "Number of harmonics specified by nh"
    "Amplitude variance explained by each harmonic is also computed"
    "returns list of amplitude variance explained, phase angles and amplitudes for each harmonic"
    pi=mh.pi
    cos=mh.cos
    sin=mh.sin
    atan=mh.atan
    mean=np.mean(obs)
    nobs=len(obs)
    S=np.var(obs)

    harm_funcs=np.arange(1,nh+1)
    #this section works out the amplitude and phaseshift
    #1st work out A1 and B1 using formulae
    list_ps=[]
    list_C=[]
    for h in harm_funcs:
        list_A=[]
        list_B=[]
        for t,y in zip(time,obs):
            A=y*(cos((2*pi*h*t)/nobs))
            B=y*(sin((2*pi*h*t)/nobs))
            list_A.append(A)
            list_B.append(B)
        A_total=np.sum(list_A)
        B_total=np.sum(list_B)
        Ax=(2*(A_total))/nobs
        Bx=(2*(B_total))/nobs

        #defining C as the amplitude
        C=np.sqrt((Ax**2 + Bx**2))

        if Ax>0:
            phi=atan(Bx/Ax)
        elif Ax<0:
            if 0 < atan(Bx/Ax)+pi < 2*pi:
                phi=atan(Bx/Ax)+pi
            else:
                phi=atan(Bx/Ax)-pi
        elif Ax==0:
            phi=pi/2

        #convert phase angle in radians to degrees.
        ps_deg=mh.degrees(phi)
        list_ps.append(ps_deg)
        list_C.append(C)
    sum_amps=np.sum(list_C)
    ex_var_list=[]
    amp_var_list=[]


    for h in harm_funcs:
        index=h-1
        x=list_C[index]
        amp_var=(x/sum_amps)*100
        amp_var_list.append(amp_var)
        ex_variance=(((nobs/2)*(x**2))/((nobs-1)*S))*100

        ex_var_list.append(ex_variance)

    return list_C, list_ps, ex_var_list, amp_var_list