;
;
;
;
;


n_step = 3L

t_c = fltarr(n_step)
t_g = fltarr(n_step)
par = fltarr(n_step)



FOR i=0L, n_step-1 DO BEGIN

   np = 512L * 2L^i

   par[i] = np

   ns = 4L^(n_step - i -1)

   print, '>>> Computing for ', np, ' particles...', ns

   ;continue

   x    = randomu(1234567, np)
   y    = randomu(2345678, np)
   z    = randomu(3456789, np)
   mass = 1.0

   t0 = systime(1)
   print, '>>> CUDA', t0
   FOR j=0, ns-1 DO up_g = potential_nbody_gpu(x,y,z,1.0)
   t1 = systime(1)
   t_g[i] = (t1-t0)/float(ns)
   print, '   t1 ',t1

   t0 =systime(1)
   print, '>>> C...',t0
   FOR j=0, ns-1 DO up_c = halo_potential_c(x,y,z,1.0)
   t1 =systime(1)
   t_c[i] = (t1-t0)/float(ns)
   print, '   t1 ', t1

   ;print, '>>> IDL...'
   ;up_i = halo_potential_slow(x3,y3,z3) ; * 8.5532206e+37

   
ENDFOR

stop

print, par
print, t_c/t_g

set_plot, 'ps'
device, /encapsulated, xsize=10, ysize=8, filename='potentital_c_v_cuda.eps'
plot, par, t_c/t_g, psym=10, xthick=2,ythick=2, thick=2, xtitle='!6N!Dp!N', ytitle='t!DC!N / t!DCUDA!N',$
      charsize=1.0, title = 'Direct potential: C vs CUDA',/xlog,xrange=[500,20000],xs=1
device,/close
set_plot, 'x'


end
