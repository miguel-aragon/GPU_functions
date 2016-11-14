;
;
;
;  INPUT:
;    x,y,z     [kpc]
;    mass      [kpc]
;    softening [kpc]
;
;  OUTPUT:
;    U_P    potential energy array
;
;
;
;-
FUNCTION potential_nbody_ccc, x,y,z, mass, softening

    ;--- Number of particles (C_INT)
    n_part = long(N_ELEMENTS(x))

    ;--- Make sure we pass the right type (C_DOUBLE)
    x    = float(x)
    y    = float(y)
    z    = float(z)

    ;--- Declare potential energy array. IDL must do the allocation
    ;    for all input and output arrays passed between IDL-C
    U_p = fltarr(n_part)

    ;--- Define external function name
    ext_function = 'potential_nbody_ccc'

    ;--- Cal function
    result = CALL_EXTERNAL(ext_function+'.so', $
			   ext_function, x,y,z, n_part, softening, U_p, $
			   VALUE=[0,0,0,1,1,0], /F_VALUE, /CDECL, /AUTO_GLUE)


    ;c_mass = double(1.989E33)   ;--- [M_o -> gm]
    ;c_dist = double(3.086E21)   ;--- [Kpc -> cm]
    ;G_grav = double(6.672E-8)   ;--- [ cm^3 g^(-1) s^(-2) ]
    ;C_norm = 8.5532206e+37 ;--- c_mass^2 * G_grav / c_dist

    return, U_p

END


