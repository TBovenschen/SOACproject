#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 10:52:56 2020

@author: tychobovenschen
"""

function [z0,ftauS,d,PsiH] = drag_model(H,lambda,slope,model,model_Cr)
# Calculates the roughness length for momentum of a surface, given the height and
# spacing density of the roughness obstacles. Different bulk drag models can be used.
#
#   Input:
#   - H:        height of obstacles (m)
#   - lambda :  spacing density of obstacles (-)
#   - slope:    mean slope of obstacles (-) NOT USED
#   - model:    type of model. Can be 'Lettau1969', 'Raupach1992', 'Raupach1994', 'Macdonald1998'
#   - model_Cr: parametrization of form drag coefficient. Can be
#   'constant', 'Banke1980', 'Garbrecht2002', 'KeanSmith2006' 
#
#   Output:
#   - z0:       aerodynamic roughness length (for momentum) (m)
#   - ftauS:    fraction of skin friction of total drag (-)
#   - d:        displacement height (m)
#   - PsiH:     correction of horizontal wind profile at z=H (-)
#
# Example:
#   [z0,ftauS,d,PsiH] =
#   drag_model(1.1,0.13,0,'Raupach1992','Garbrecht2002')
#
# References:
# ï»¿- Lettau H (1969)
#   Note on Aerodynamic Roughness-Parameter Estimation on the Basis of Roughness-Element Description.
#   J. Appl. Meteor. 8:828â832
# ï»¿- Raupach MR (1992) 
#   Drag and drag partition on rough surfaces.
#   Boundary-Layer Meteorol 60:375â395
# ï»¿- Raupach MR (1994) 
#   Simplified expressions for vegetation roughness length and zero-plane displacement as functions of canopy height and area index. 
#   Boundary-Layer Meteorol 71:211â216. 
# ï»¿- Macdonald RW, Griffiths RF, Hall DJ (1998) 
#   An improved method for the estimation of surface roughness of obstacle arrays. 
#   Atmos Environ 32:1857â1864.
#
# history:
#   - 29/09/2020: written by Maurice van Tiggelen (UU) for the SOAC project 


# Global constants
kappa   = 0.4 #von karman constant

# Parameters
Cs10    = 0.00083 # drag coefficient for skin friction at z=10m (default = 0.00083)
c       = 0.25 # parameter of Raupach(1992) model
c1      = 7.5 # parameter of Raupach(1992) model
A       = 4 # parameter of Macdonald(1998) model

# Initialize output
z0      = nan
ftauS   = nan
d       = nan
PsiH    = nan

# Parametrization of form drag coefficient (not used in Lettau1969)
switch model_Cr
    case 'constant'
        Cr = 0.1
    case 'Banke1980'
        Cr = (0.012 + 0.012.*slope)
    case 'Garbrecht2002'
        if H<2.5527
            Cr = (0.185 + (0.147*H))/2
        else
            Cr = 0.11*log(H/0.2)
        end
    case 'KeanSmith2006'
        Cr = 0.8950*exp(-0.77*(0.5./(4*lambda)))
end

# Drag model
switch (model):
    case 'Lettau1969' 
        # roughness length
        z0 = 0.5.*H*lambda
                    
    case 'Raupach1992'
        # displacement height
        d = H*(1 - (1-exp(-(c1.*lambda).^(0.5)))/(c1.*lambda).^0.5)
        
        # Profile correction
        PsiH = 0.193

        # Drag coefficient for skin friction at z=H
        Cs = (Cs10.^(-0.5) - (1/kappa).*(log((10-d)/(H-d))-PsiH)).^(-2)
        
        a = (c*lambda/2)*(Cs+lambda*Cr).^(-0.5)
        
        # model for U/u*
        X = a
        crit=1e5
        while(crit>1e-12):
            Xold=X
            X = a*exp(X)
            crit=abs((X-Xold)./(X))
        gamma = 2*X./(c*lambda)

        # roughness length
        z0 = (H-d).*exp(-kappa.*gamma).*exp(PsiH)
        
        # stress partionning
        beta  = Cr./Cs
        ftauS = 1./(1+beta*lambda)
                
    case 'Raupach1994'
        # displacement height
        d = H*(1 - (1-exp(-(c1.*lambda).^(0.5)))/(c1.*lambda).^0.5)
        
        PsiH = 0.193
        
        # Drag coefficient for skin friction at z=H
        Cs = (Cs10.^(-0.5) - (1/kappa).*(log((10-d)/(H-d))-PsiH)).^(-2)
        
        # model for U/u*
        gamma_s = (Cs + Cr.*lambda).^(-0.5)
        
        # roughness length
        z0 = (H-d).*(exp(kappa.*gamma_s)-PsiH).^(-1)

        # stress partionning
        beta  = Cr./Cs
        ftauS = 1./(1+beta*lambda)
        
    case 'Macdonald1998'
        # displacement height
        d =  1 + A.^(-lambda) .* (lambda  - 1)
        # roughness length
        z0 = (H - d).*exp(-((Cr./(kappa.^2)).*lambda) .^(-0.5))   
end

end

