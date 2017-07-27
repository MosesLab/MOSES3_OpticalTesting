;Test surface_height2SLOPE by against a circular mirror surface_height with known RMS slope.
;
;To run the test:
;  IDL>  .run surface_height2slope_test
;
;CCK 2016-Jun-20

;First, produce a 'perfect' mirror (surface_height=0).
N = 512               ; NxN = square array size
mirror_diameter = 200 ; mm
dx = 1.1*mirror_diameter/N ; mm --- This will be the pixel size
lambda = 650e-6        ; wavelength, mm
mirror_radius_pixels = mirror_diameter/(2.0*dx)
surface_height = replicate(0.0,N,N)
r_coordinate = shift(dist(N), N/2, N/2)
off_mirror = where( r_coordinate gt mirror_radius_pixels)
NaN = -0.0/0.0        ; One way of coding IEEE NaN
surface_height(off_mirror) =  NaN;

;Now add a ripple of known frequency.
x = (findgen(N) * dx) # replicate(1.0,N) ;mirror surface_height x-coordinate
amplitude = 0.033        ;waves
period = mirror_diameter/10.0 ;ripple period, mm
k = 2 * !pi / period   ;ripple wavenumber, radians per mm
ripple = amplitude * sin(k*x) ;ripple in waves
ripple_RMS = lambda * amplitude / sqrt(2.0) ;analytic RMS ripple, mm
surface_height += ripple          ;apply ripple to surface_height
slope = lambda * amplitude * k * cos(k*x) ;analytic slope in RADIANS
slope_RMS = lambda * amplitude * k / sqrt(2.0) ;analytic RMS slope, radians
print,'Analytic RMS slope (radians):   ', slope_RMS
print,'Confirm analytic result:        ', sqrt(mean(slope^2))

;Display our test image
window, 1, xsize=N, ysize=N, title='surface_height2SLOPE Test Image'
loadct,0
tv, displayscale(surface_height, /ct)
loadct,0

;See if surface_height2slope does the job!
integration_length = period/10.0 ;mm
filter_length = integration_length/2.0
slope_RMS_numerical = wfe2slope(surface_height, integration_length, filter_length, $
   dx=dx, lambda=lambda)
print,'surface_height2slope() RMS slope (radians):', slope_RMS_numerical
print

;Try varying the filter_length and integration_length
mx = 0.25*10^(findgen(50)/20)
   ;multipliers for existing filter & integration length parameters
Ntrials = n_elements(mx)
integration_length_arr = mx * integration_length
filter_length_arr = mx * filter_length
slope_RMS_arr = replicate(slope_RMS, Ntrials)
slope_RMS_num_arr = findgen(Ntrials) ;storage for results
for i=0, Ntrials-1 do begin
   slope_RMS_num_arr[i] = wfe2slope(surface_height, integration_length_arr[i], $
      filter_length_arr[i], dx=dx, lambda=lambda)
endfor
window, 2, title='slope error vs measurement prams'
plot, integration_length_arr, 1e6*slope_RMS_num_arr, /xlog, charsize=2, $
   ytitle = 'RMS slope (urad)', $
   xtitle = 'integration length (mm) = 2 x filter length (mm)'
oplot, integration_length_arr, 1e6*replicate(slope_RMS, Ntrials), linestyle=2
oplot, [1,1]*period/2, [0.001,1000], linestyle=1 ;plot semi-period
oplot, [1,1]*dx, [1e-3,1000], linestyle=3 ;plot pixel size
al_legend, ['calculated with surface_height2slope()','exact','1/2 ripple period', $
   'pixel size'], $
   linestyle=[0,2,1,3], charsize=2, /bottom, /left

;Try varying the angle of the ripple pattern.
thetas = 10.0*findgen(10) ;angles in degrees
Nthetas = n_elements(thetas)
slope_RMS_arr2 = replicate(slope_RMS, Nthetas)
slope_RMS_num_arr2 = findgen(Nthetas)
for i=0, Nthetas-1 do begin
   surface_height2 = rot(surface_height, thetas[i], missing=NaN, cubic=-0.5)
   slope_RMS_num_arr2[i] = wfe2slope(surface_height2, integration_length, filter_length, $
      dx=dx, lambda=lambda)
endfor
window, 3, title='slope error vs clocking of ripple pattern'
plot, thetas, 1e6*slope_RMS_num_arr2, charsize=2, $
   xtitle='ripple angle (degrees clockwise)', $
   ytitle = 'RMS slope (urad)'
oplot, thetas, 1e6*replicate(slope_RMS, Ntrials), linestyle=2
al_legend, ['calculated with surface_height2slope()','exact'], $
   linestyle=[0,2], charsize=2, /bottom, /left


end