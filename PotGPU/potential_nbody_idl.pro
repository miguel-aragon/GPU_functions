function potential_nbody_idl, x_in,y_in,z_in,  mass, softening
;
;
;  INPUT:
;    xyz  positions [kpc]
;    mass           [M_o] 
;    softening      [kpc]
;

;c_mass = double(1.989E33)  ;--- [M_o -> gm]
;c_dist = double(3.086E21)  ;--- [Kpc -> cm]
;G_grav = double(6.672E-8)  ;--- [ cm^3 g^(-1) s^(-2) ]

;--- Particles in this halo
n_sub = N_ELEMENTS(x_in)

;--- Convert units
x = float(x_in)
y = float(y_in)
z = float(z_in)

;--- squared softening
soft2 = softening^2

;--- Array of potential energy
U_p = fltarr(n_sub)

;--- Compute potential energy the hard way...
FOR i=0L, n_sub-1 DO BEGIN
    ;--- Avoid itself
    sub = where(x NE x[i])
    ;--- Full potential
    U_p[i] = total(1.0/sqrt((x[sub]-x[i])^2 + (y[sub]-y[i])^2 + (z[sub]-z[i])^2 + soft2 ) )
ENDFOR

return, U_p

end

