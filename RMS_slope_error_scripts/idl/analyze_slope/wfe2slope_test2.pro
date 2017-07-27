;Test WFE2SLOPE by against a circular mirror WFE with known RMS slope
;NOTE THAT I HAVE TESTED WFE SLOPE, NOT MIRROR SLOPE (WHICH WOULD BE SMALLER /2).
;
;Also tests RMS ripple calculation using rmsripple().
;
;CCK 2015-Mar-11

;First, produce a 'perfect' mirror (wfe=0).
mirror_diameter = 200 ; mm
lambda = 30e-6        ; wavelength, mm
N = 512               ; NxN = square array size
wfe = mirror(mirror_diameter, 1e-30, dx=dx, N=N, /circular) 
   ;Note that mirror() defaults to nonzero wfe when rmswfe=0.

;Now add a ripple of known frequency.
x = (findgen(N) * dx) # replicate(1.0,N) ;mirror surface x-coordinate
amplitude = 1.0        ;waves
period = mirror_diameter/20.0 ;ripple period, mm
k = 2 * !pi / period   ;radians per mm
ripple = amplitude * sin(k*x) ;ripple in waves
ripple_RMS = lambda * amplitude / sqrt(2.0) ;analytic RMS ripple, mm
wfe += ripple          ;apply ripple to WFE
slope = lambda * amplitude * k * cos(k*x) ;analytic slope in RADIANS
slope_RMS = lambda * amplitude * k / sqrt(2.0) ;analytic RMS slope, radians
print,'Analytic RMS slope (radians):   ', slope_RMS
print,'Confirm analytic result:        ', sqrt(mean(slope^2))

;See if wfe2slope does the job!
integration_length = period/10.0 ;mm
filter_length = integration_length/2.0
slope_RMS_numerical = wfe2slope(wfe, integration_length, filter_length, $
   dx=dx, lambda=lambda)
print,'wfe2slope() RMS slope (radians):', slope_RMS_numerical
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
   slope_RMS_num_arr[i] = wfe2slope(wfe, integration_length_arr[i], $
      filter_length_arr[i], dx=dx, lambda=lambda)
endfor
window, 1, title='slope error vs measurement prams'
plot, integration_length_arr, 1e6*slope_RMS_num_arr, /xlog, charsize=2, $
   ytitle = 'RMS slope (urad)', $
   xtitle = 'integration length (mm) = 2 x filter length (mm)'
oplot, integration_length_arr, 1e6*replicate(slope_RMS, Ntrials), linestyle=2
oplot, [1,1]*period/2, [0.001,1000], linestyle=1 ;plot semi-period
oplot, [1,1]*dx, [3,1000], linestyle=3 ;plot pixel size
al_legend, ['calculated with wfe2slope()','exact','1/2 ripple period', $
   'pixel size'], $
   linestyle=[0,2,1,3], charsize=2, /bottom, /left

;Try varying the angle of the ripple pattern.
thetas = 10.0*findgen(10) ;angles in degrees
Nthetas = n_elements(thetas)
slope_RMS_arr2 = replicate(slope_RMS, Nthetas)
slope_RMS_num_arr2 = findgen(Nthetas)
for i=0, Nthetas-1 do begin
   wfe2_real = rot(float(wfe), thetas[i], missing=0.0, cubic=-0.5)
   wfe2 = complex(wfe2_real, imaginary(wfe))
   slope_RMS_num_arr2[i] = wfe2slope(wfe2, integration_length, filter_length, $
      dx=dx, lambda=lambda)
endfor
window, 2, title='slope error vs clocking of ripple pattern'
plot, thetas, 1e6*slope_RMS_num_arr2, charsize=2, $
   xtitle='ripple angle (degrees clockwise)', $
   ytitle = 'RMS slope (urad)'
oplot, thetas, 1e6*replicate(slope_RMS, Ntrials), linestyle=2
al_legend, ['calculated with wfe2slope()','exact'], $
   linestyle=[0,2], charsize=2, /bottom, /left

;Now test our numerical method for ripple:
;N.B.: I am looking at the wavefront (not surface) ripple.
ripple_longest_period = period*3.0
ripple_bandwidth = 10.0
ripple_RMS_numerical = rmsripple(wfe, ripple_longest_period, ripple_bandwidth, $
   dx=dx, lambda=lambda)
print,'Analytic RMS ripple (nm):', 1e6 * ripple_RMS
print,'Confirm analytic result: ', 1e6 * lambda * sqrt(mean(ripple^2))
print,'rmsripple() result:      ', 1e6 * ripple_RMS_numerical
end