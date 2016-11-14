;
;
;
;
;



   np = 2010L

   x    = randomu(1234567, np)
   y    = randomu(2345678, np)
   z    = randomu(3456789, np)
   mass = 1.0

   ;softening = 0.0001

   print, '>>> CUDA'
   up_g = potential_nbody_gpu(x,y,z,1.0, softening)


   print, '>>> C...'
   up_c = potential_nbody_ccc(x,y,z,1.0, softening)


   print, '>>> IDL...'
   up_i = potential_nbody_idl(x,y,z,1.0, softening)



 ; * 8.5532206e+37

   print, up_i[0:5]


   print, '---------------------'
   more, up_g
   more, up_c
   more, up_i
   print, '---------'
   more, up_g - up_c
   more, up_g - up_i
   print, '---------'

   print, 't(g-c)      = ', total(abs(up_g - up_c))
   print, 't(g-i)      = ', total(abs(up_g - up_i))
   
   print, '---------'
   print, 't(g) - t(c) = ', total(up_g) - total(up_c)
   print, 't(g) - t(i) = ', total(up_g) - total(up_i)



end
