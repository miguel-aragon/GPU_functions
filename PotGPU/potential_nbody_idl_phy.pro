function potential_nbody_idl_phy, x_in,y_in,z_in,  mass, softening
;
;
;  INPUT:
;    xyz  positions [kpc]
;    mass           [M_o] 
;    softening      [kpc]
;

c_mass = double(1.989E33)  ;--- [M_o -> gm]
c_dist = double(3.086E21)  ;--- [Kpc -> cm]
G_grav = double(6.672E-8)  ;--- [ cm^3 g^(-1) s^(-2) ]

;--- Particles in this halo
n_sub = N_ELEMENTS(x_in)

;--- Convert units
x = double(x_in)*c_dist
y = double(y_in)*c_dist
z = double(z_in)*c_dist

;--- squared softening
soft2 = softening^2

;--- Array of potential energy
U_p = dblarr(n_sub)

;--- Compute potential energy the hard way...
FOR i=0L, n_sub-1 DO BEGIN
    ;--- Avoid itself
    sub = where(x NE x[i])
    ;--- Full potential
    U_p[i] = G_grav * c_mass^2 * total(1.0/sqrt((x[sub]-x[i])^2 + (y[sub]-y[i])^2 + (z[sub]-z[i])^2 + soft2 ) )
ENDFOR

return, U_p

end

