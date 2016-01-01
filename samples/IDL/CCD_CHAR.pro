PRO CCD_CHAR

SUM = dialog_pickfile(Title="Select SUM Images", /multiple_files)     ; Reading sum, difference, and
DIF = dialog_pickfile(Title="Select DIF Images", /multiple_files)     ; bias files in. These should
BIAS = dialog_pickfile(Title="Select BIAS Images", /multiple_files)   ; have been created beforehand.

PRINT, "Processing SUM images ..."
SUMpix = DBLARR(size(SUM, /N_ELEMENTS))
SUM_mean = DBLARR(size(SUM, /N_ELEMENTS))
SUM_stddev = DBLARR(size(SUM, /N_ELEMENTS))
SUM_sigmean = DBLARR(size(SUM, /N_ELEMENTS))

for i=0,((size(SUM, /N_ELEMENTS))-1) DO BEGIN
  PRINT, "Processing SUM image ",strtrim(i+1,2)," of ",strtrim((size(SUM, /N_ELEMENTS)),2),"."
  image_display,SUM[i],area,mean,stddev,sig_mean
  SUMpix[i] = area                                  ; Generating statistics of the sum images.
  SUM_mean[i] = mean
  SUM_stddev[i] = stddev
  SUM_sigmean[i] = sig_mean
endfor

PRINT, "SUM images completed!"

PRINT, "Processing DIF images..."
DIFpix = DBLARR(size(DIF, /N_ELEMENTS))
DIF_mean = DBLARR(size(DIF, /N_ELEMENTS))
DIF_stddev = DBLARR(size(DIF, /N_ELEMENTS))
DIF_sigmean = DBLARR(size(DIF, /N_ELEMENTS))

for i=0,((size(DIF, /N_ELEMENTS))-1) DO BEGIN
  PRINT, "Processing DIF image ",strtrim(i+1,2)," of ",strtrim((size(DIF, /N_ELEMENTS)),2),"."
  image_display,DIF[i],area,mean,stddev,sig_mean
  DIFpix[i] = area                                ; Generating statistics of the difference images.
  DIF_mean[i] = mean
  DIF_stddev[i] = stddev
  DIF_sigmean[i] = sig_mean
endfor

PRINT, "DIF images completed!"

PRINT, "Processing BIAS images..."
Xpix = DBLARR(size(BIAS, /N_ELEMENTS))
X_mean = DBLARR(size(BIAS, /N_ELEMENTS))
X_stddev = DBLARR(size(BIAS, /N_ELEMENTS))
X_sigmean = DBLARR(size(BIAS, /N_ELEMENTS))

for i=0,((size(BIAS, /N_ELEMENTS))-1) DO BEGIN
  PRINT, "Processing BIAS image ",strtrim(i+1,2)," of ",strtrim((size(BIAS, /N_ELEMENTS)),2),"."
  image_display,BIAS[i],area,mean,stddev,sig_mean
  Xpix[i] = area
  X_mean[i] = mean
  X_stddev[i] = stddev
  X_sigmean[i] = sig_mean
endfor

PRINT, "BIAS images completed!"

t_exp = [0.25, 0.5, INDGEN(30, START=1)]

PRINT, "Processing preliminary plot..."
X_meanf = DBLARR((size(SUM_mean, /N_ELEMENTS)))
for i=0,((size(X_mean, /N_ELEMENTS))-1) DO X_meanf[i] = X_mean[i]

x = SUM_mean - 2*X_meanf
y = DIF_stddev^2
xtitle1 = "Measured Counts (counts/px)"
ytitle1 = "CCD Noise (counts^2/px^2)"

xerr_i = sqrt(SUM_sigmean^2+4*(mean(X_sigmean)^2))
xerr = DBLARR((size(x, /N_ELEMENTS)))
for i=0,((size(xerr_i, /N_ELEMENTS))-1) DO xerr[i] = xerr_i[i]

SET_PLOT, 'PS'
DEVICE, FILE='measure1.eps', /ENCAPSULATED
plot, t_exp, y, psym=7, xtitle="Exposure Time (s)", ytitle=ytitle1, yrange=[0,3.5d4];, $
DEVICE, /CLOSE
DEVICE, FILE='measure2.eps', /ENCAPSULATED
ploterror, x, y, xerr, psym=7, xtitle=xtitle1, ytitle=ytitle1, yrange=[0,3.5d4];, $
DEVICE, /CLOSE

PRINT, "Fitting linear model..."

dat1 = 14
fit1 = LINFIT(x[0:dat1], y[0:dat1], SIGMA=sig_lin1, CHISQR=chisqr1, PROB=prob1)
model1 = fit1[0] + fit1[1]*x
Err_i = (x^2)*(sig_lin1[1]^2)+(fit1[1]^2)(xerr^2)
Err1 = (sqrt(sig_lin1[0]^2 + Err_i))

dat2 = 13
fit2 = LINFIT(x[0:dat2], y[0:dat2], SIGMA=sig_lin2, CHISQR=chisqr2, PROB=prob2)
model2 = fit2[0] + fit2[1]*x
Err_i = (x^2)*(sig_lin2[1]^2)+(fit2[1]^2)(xerr^2)
Err2 = (sqrt(sig_lin2[0]^2 + Err_i))

dat3 = 15
fit3 = LINFIT(x[0:dat3], y[0:dat3], SIGMA=sig_lin3, CHISQR=chisqr3, PROB=prob3)
model3 = fit3[0] + fit3[1]*x
Err_i = (x^2)*(sig_lin3[1]^2)+(fit3[1]^2)(xerr^2)
Err3 = (sqrt(sig_lin3[0]^2 + Err_i))

DEVICE, FILE='linfit.eps', /ENCAPSULATED
ploterror, x, y, xerr, psym=7, xtitle=xtitle1, ytitle=ytitle1, yrange=[0,5.0d4];, $
oploterror, x, model1, Err1, linestyle=0
DEVICE, /CLOSE
DEVICE, FILE='linfit_multi.eps', /ENCAPSULATED
ploterror, x, y, xerr, psym=7, xtitle=xtitle1, ytitle=ytitle1, yrange=[0,5.0d4];, $
oploterror, x, model1, Err1, linestyle=0
oploterror, x, model2, Err2, linestyle=2
oploterror, x, model3, Err3, linestyle=3
DEVICE, /CLOSE

deviation1 = model1 - y
deviation2 = model2 - y
deviation3 = model3 - y

zero = INTARR(size(x, /N_ELEMENTS))
DEVICE, FILE='lin_dev.eps', /ENCAPSULATED
plot, t_exp, deviation1, xtitle="Exposure Time (s)", ytitle=ytitle1, linestyle=0;, $
oplot, t_exp, deviation2, linestyle=2
oplot, t_exp, deviation3, linestyle=3
oplot, t_exp, zero, thick=2
DEVICE, /CLOSE

DEVICE, FILE='lin_dev_f.eps', /ENCAPSULATED
plot, t_exp, deviation1, xtitle="Exposure Time (s)", psym=7, ytitle=ytitle1, linestyle=2
oplot, t_exp, zero, thick=2
DEVICE, /CLOSE

PRINT, "Characterising CCD..."

GAIN = 1.0d0/fit1[1]
GAINerr = (1.0d0/fit1[1]^2)*(sig_lin1[1])
RN = sqrt(fit1[0]/2)*GAIN
sig1 = ((GAIN^2*sig_lin1[0]^2)/(8*fit1[0]))+(fit1[0]*GAINerr)/2
RNerr = sqrt(sig1)
FWD = (x[dat1]/2)*GAIN
FWDerr = sqrt((GAIN^2)*((xerr[dat1]/2)^2)+((x[dat1]/2)^2)*(GAINerr^2))

PRINT, "For FIT1 with ",strtrim(dat1,2)," data points..."
PRINT, "Gain: ",strtrim(GAIN,2)," electrons/ADU, Error: ",strtrim(GAINerr,2),"."
PRINT, "Read Noise: ",strtrim(RN,2)," electrons/px, Error: ",strtrim(RNerr,2),"."
PRINT, "Pixel Full-Well Depth: ",strtrim(FWD,2)," electrons/px, Error: ",strtrim(FWDerr,2),"."

GAIN = 1.0d0/fit2[1]
GAINerr = (1.0d0/fit2[1]^2)*(sig_lin2[1])
RN = sqrt(fit2[0]/2)*GAIN
sig1 = ((GAIN^2*sig_lin2[0]^2)/(8*fit2[0]))+(fit2[0]*GAINerr)/2
RNerr = sqrt(sig1)
FWD = (x[dat2]/2)*GAIN
FWDerr = sqrt((GAIN^2)*((xerr[dat2]/2)^2)+((x[dat2]/2)^2)*(GAINerr^2))

PRINT, "For FIT2 with ",strtrim(dat2,2)," data points..."
PRINT, "Gain: ",strtrim(GAIN,2)," electrons/ADU, Error: ",strtrim(GAINerr,2),"."
PRINT, "Read Noise: ",strtrim(RN,2)," electrons/px, Error: ",strtrim(RNerr,2),"."
PRINT, "Pixel Full-Well Depth: ",strtrim(FWD,2)," electrons/px, Error: ",strtrim(FWDerr,2),"."

GAIN = 1.0d0/fit3[1]
GAINerr = (1.0d0/fit3[1]^2)*(sig_lin3[1])
RN = sqrt(fit3[0]/2)*GAIN
sig1 = ((GAIN^2*sig_lin3[0]^2)/(8*fit3[0]))+(fit3[0]*GAINerr)/2
RNerr = sqrt(sig1)
FWD = (x[dat3]/2)*GAIN
FWDerr = sqrt((GAIN^2)*((xerr[dat3]/2)^2)+((x[dat3]/2)^2)*(GAINerr^2))

PRINT, "For FIT3 with ",strtrim(dat3,2)," data points..."
PRINT, "Gain: ",strtrim(GAIN,2)," electrons/ADU, Error: ",strtrim(GAINerr,2),"."
PRINT, "Read Noise: ",strtrim(RN,2)," electrons/px, Error: ",strtrim(RNerr,2),"."
PRINT, "Pixel Full-Well Depth: ",strtrim(FWD,2)," electrons/px, Error: ",strtrim(FWDerr,2),"."

end
