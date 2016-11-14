;
;
;
;  INPUT:
;    x,y,z     positions [kpc]
;    mass                [M_o]
;    softening           [kpc]
;
;  OUTPUT:
;    U_P    potential energy array
;
;
;
;-
FUNCTION potential_nbody_gpu, x,y,z, mass, softening

    ;--- library file
    lib_file = 'potential_nbody_gpu.so'

    ;--- Number of particles (C_INT)
    n_part = long(N_ELEMENTS(x))

    ;--- Make sure we pass the right type (C_FLOAT)
    pos      = fltarr(3,n_part)
    pos[0,*] = float(x)
    pos[1,*] = float(y)
    pos[2,*] = float(z)
    mass     = fltarr(n_part)+1.0

    ;--- Set softening
    dump = CALL_EXTERNAL(lib_file, $
                        'set_softening', float(softening), $
                        VALUE=[1], /F_VALUE, /CDECL, /AUTO_GLUE)


    ;--- Init computation
    dump = CALL_EXTERNAL(lib_file, $
			'init_nbody', n_part, pos, mass, $
                        VALUE=[1,0,0], /F_VALUE, /CDECL, /AUTO_GLUE)

    ;--- Compute potential
    dt = float(1.0)
    dump = CALL_EXTERNAL(lib_file, $
			'compute_potential',  $
                        /F_VALUE, /CDECL, /AUTO_GLUE)


    ;--- Retrieve data
    Poten = fltarr(n_part)
    dump = CALL_EXTERNAL(lib_file, $
			'get_potential', Poten, $
                        VALUE=[0], /F_VALUE, /CDECL, /AUTO_GLUE)

    ;--- Clean GPU memory
    dump = CALL_EXTERNAL(lib_file, $
			'dealloc_nbody', $
                        /F_VALUE, /CDECL, /AUTO_GLUE)



    ;c_mass = double(1.989E33)   ;--- [M_o -> gm]
    ;c_dist = double(3.086E21)   ;--- [Kpc -> cm]
    ;G_grav = double(6.672E-8)   ;--- [ cm^3 g^(-1) s^(-2) ]
    ;C_norm = 8.5532206e+37 ;--- c_mass^2 * G_grav / c_dist

    return, Poten

END


