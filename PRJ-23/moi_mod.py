'''
Script to approximate moments of inertia
Ney Rafael Sêcco
ney@ita.br
Instituto Tecnológico de Aeronáutica
09-2022
'''

# IMPORTS
import numpy as np
import matplotlib.pyplot as plt

def surf_moi(xr,yr,zr,cr,tcr,
             xt,yt,zt,ct,tct,
             xref,yref,zref,
             mass,
             Nc=100,Nb=200):
    '''
    This function computes the moment of inertia
    of a half-surface.

    Returns
    -------
    None.

    '''
    
    # Determine leading and trailing edge points
    xrLE = xr
    yrLE = yr
    zrLE = zr
    xrTE = xr+cr
    yrTE = yr
    zrTE = zr
    
    xtLE = xt
    ytLE = yt
    ztLE = zt
    xtTE = xt+ct
    ytTE = yt
    ztTE = zt
    
    # Inverse matrix for bilinear interpolation
    Ainv = np.linalg.inv(np.array([[1,0,0,0],
                                   [1,1,0,0],
                                   [1,0,1,0],
                                   [1,1,1,1]]))
    
    # Create a mesh for properties varying along the wing.
    # We do this by creating a linear interpolation of properties.
    # uu=0 at the LE and uu=1 at TE
    # vv=0 at the root and vv=1 at the tip
    uu = np.linspace(0,1,Nc)
    vv = np.linspace(0,1,Nb)
    
    UU, VV = np.meshgrid(uu, vv)
    
    # X interpolation
    iw = Ainv.dot([xrLE,xrTE,xtLE,xtTE])
    XX = iw[0] + iw[1]*UU + iw[2]*VV + iw[3]*UU*VV
    
    # Y interpolation
    iw = Ainv.dot([yrLE,yrTE,ytLE,ytTE])
    YY = iw[0] + iw[1]*UU + iw[2]*VV + iw[3]*UU*VV
    
    # Z interpolation
    iw = Ainv.dot([zrLE,zrTE,ztLE,ztTE])
    ZZ = iw[0] + iw[1]*UU + iw[2]*VV + iw[3]*UU*VV
    
    # absolute thickness interpolation
    iw = Ainv.dot([tcr*cr,tcr*cr,tct*ct,tct*ct])
    TT = iw[0] + iw[1]*UU + iw[2]*VV + iw[3]*UU*VV
    
    # Compute areas of each quadrilateral using cross product
    # of diagonals
    dX1 = XX[1:,1:]  - XX[:-1,:-1]
    dX2 = XX[1:,:-1] - XX[:-1,1:]
    dY1 = YY[1:,1:]  - YY[:-1,:-1]
    dY2 = YY[1:,:-1] - YY[:-1,1:]
    dZ1 = ZZ[1:,1:]  - ZZ[:-1,:-1]
    dZ2 = ZZ[1:,:-1] - ZZ[:-1,1:]
    
    # Cross product and area
    dXv = dY1*dZ2 - dY2*dZ1
    dYv = dZ1*dX2 - dX1*dZ2
    dZv = dX1*dY2 - dY1*dX2
    
    AA = 0.5*np.sqrt(dXv**2 + dYv**2 + dZv**2)
    
    # Get average thickness and centroid of each quadrilateral
    Xc = 0.25*(XX[:-1,:-1] + XX[1:,:-1] + XX[1:,1:] + XX[:-1,1:])
    Yc = 0.25*(YY[:-1,:-1] + YY[1:,:-1] + YY[1:,1:] + YY[:-1,1:])
    Zc = 0.25*(ZZ[:-1,:-1] + ZZ[1:,:-1] + ZZ[1:,1:] + ZZ[:-1,1:])
    Tc = 0.25*(TT[:-1,:-1] + TT[1:,:-1] + TT[1:,1:] + TT[:-1,1:])
    
    # Compute volume of each element and of the wing
    Vc = AA*Tc
    vol = np.sum(Vc)
    
    # Wing density
    rho = mass/vol
    
    # Moments of inertia
    Ixx = np.sum(rho*Vc*((Yc-yref)**2 + (Zc-zref)**2))
    Iyy = np.sum(rho*Vc*((Xc-xref)**2 + (Zc-zref)**2))
    Izz = np.sum(rho*Vc*((Xc-xref)**2 + (Yc-yref)**2))
    Ixy = np.sum(rho*Vc*(Xc-xref)*(Yc-yref))
    Ixz = np.sum(rho*Vc*(Xc-xref)*(Zc-zref))
    Iyz = np.sum(rho*Vc*(Yc-yref)*(Zc-zref))
    
    return Ixx,Iyy,Izz,Ixy,Ixz,Iyz

#======================================

def cyl_moi(xn,yn,zn,
            Lcyl,Dcyl,
            xref,yref,zref,
            mass,
            Nr=200,Nt=200,Nl=10):
    '''
    This function computes the moment of inertia
    of a half-surface.

    Returns
    -------
    None.

    '''
    
    # Compute volume of the cylinder
    vol = np.pi*Dcyl**2/4*Lcyl
    rho = mass/vol
    
    # Compute length of each disk
    tdisk = Lcyl/Nl
    
    # First we get the distribution of one disk
    # perpendicular to the x axis
            
    # Create a mesh for properties varying along the disk.
    # We do this by creating a linear interpolation of properties.
    # uu=0 at the center and uu=1 at the rim
    # vv=0 at for theta=0 and vv=1 at for theta=2*pi, starting at the left
    uu = np.linspace(0,Dcyl/2,Nr)
    vv = np.linspace(0,2*np.pi,Nt)
    
    UU, VV = np.meshgrid(uu, vv)
    
    YY = UU*np.cos(VV) + yn
    ZZ = UU*np.sin(VV) + zn
    XX = np.zeros_like(ZZ)
    
    # Compute areas of each quadrilateral using cross product
    # of diagonals
    dX1 = XX[1:,1:]  - XX[:-1,:-1]
    dX2 = XX[1:,:-1] - XX[:-1,1:]
    dY1 = YY[1:,1:]  - YY[:-1,:-1]
    dY2 = YY[1:,:-1] - YY[:-1,1:]
    dZ1 = ZZ[1:,1:]  - ZZ[:-1,:-1]
    dZ2 = ZZ[1:,:-1] - ZZ[:-1,1:]
    
    # Cross product and area
    dXv = dY1*dZ2 - dY2*dZ1
    dYv = dZ1*dX2 - dX1*dZ2
    dZv = dX1*dY2 - dY1*dX2
    
    AA = 0.5*np.sqrt(dXv**2 + dYv**2 + dZv**2)
    
    # Get average thickness and centroid of each quadrilateral
    Xc = 0.25*(XX[:-1,:-1] + XX[1:,:-1] + XX[1:,1:] + XX[:-1,1:])
    Yc = 0.25*(YY[:-1,:-1] + YY[1:,:-1] + YY[1:,1:] + YY[:-1,1:])
    Zc = 0.25*(ZZ[:-1,:-1] + ZZ[1:,:-1] + ZZ[1:,1:] + ZZ[:-1,1:])
    
    # Compute volume of each element and of the disk
    Vc = AA*tdisk
    
    # Moments of inertia of each disk
    Ixx = 0.0
    Iyy = 0.0
    Izz = 0.0
    Ixy = 0.0
    Ixz = 0.0
    Iyz = 0.0
    for ii in range(Nl):
        
        # Ajust X coordinate
        Xc[:] = xn + tdisk/2 + ii*tdisk
        
        Ixx = Ixx + np.sum(rho*Vc*((Yc-yref)**2 + (Zc-zref)**2))
        Iyy = Iyy + np.sum(rho*Vc*((Xc-xref)**2 + (Zc-zref)**2))
        Izz = Izz + np.sum(rho*Vc*((Xc-xref)**2 + (Yc-yref)**2))
        Ixy = Ixy + np.sum(rho*Vc*(Xc-xref)*(Yc-yref))
        Ixz = Ixz + np.sum(rho*Vc*(Xc-xref)*(Zc-zref))
        Iyz = Iyz + np.sum(rho*Vc*(Yc-yref)*(Zc-zref))
    
    return Ixx,Iyy,Izz,Ixy,Ixz,Iyz

#======================================
    
def point_moi(x,y,z,
              xref,yref,zref,
              mass):
    '''
    This function computes the moment of inertia
    of a point mass.

    Returns
    -------
    None.

    '''
    
    # Moments of inertia
    Ixx = mass*((y-yref)**2 + (z-zref)**2)
    Iyy = mass*((x-xref)**2 + (z-zref)**2)
    Izz = mass*((x-xref)**2 + (y-yref)**2)
    Ixy = mass*(x-xref)*(y-yref)
    Ixz = mass*(x-xref)*(z-zref)
    Iyz = mass*(y-yref)*(z-zref)
    
    return Ixx,Iyy,Izz,Ixy,Ixz,Iyz