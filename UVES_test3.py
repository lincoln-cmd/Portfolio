def region_around_line(w, flux, cont):
    indcont = ((w > cont[0][0]) & ( w < cont[0][1])) | ((w > cont[1][0]) & (w < cont[1][1]))
    indrange = (w > cont[0][0]) & (w < cont[1][1])
    f = np.zeros((flux.shape[0], indrange.sum()))
    #print('incont: ', indcont)
    #print('indrange : ', indrange)
    #print('f : ', f)
    #print('flux: ', flux)
    #print('indrage.sum() : ', indrange.sum())
    for i in range(flux.shape[0]):
        linecoeff = np.polyfit(w[indcont], flux[i, indcont], 2)
        f[i,:] = flux[i, indrange] / np.polyval(linecoeff, w[indrange].value)
    return w[indrange], f

wcaII, fcaII = region_around_line(wavelength, flux, [[3925*u.AA, 3930*u.AA], [3938*u.AA, 3945*u.AA]])
#print(wcaII)
#print((fcaII[0]), len(fcaII[1]), len(fcaII[-1]))

# Ca2 = wavelength : 393.3682nm, line width : 2.0253nm
# Ca2(in code) = wavelength : 393.366nm, equivalent width : 20.21238214515653
ew = fcaII[0,:] - 1.
ew = ew[:-1] * np.diff(wcaII.to(u.AA).value)
#print(ew.sum())

delta_lam = np.diff(wcaII.to(u.AA).value)
ew = np.sum((fcaII - 1.)[:,:-1] * delta_lam[np.newaxis, :], axis = 1)
#print(ew)

from astropy.table import Column, Table
from astropy.io import ascii

datecol = Column(name = 'Obs Date', data = date)
pcol = Column(name = 'phase', data = delta_p, format = '{:.1f}')
ewcol = Column(name = 'EW', data = ew, format = '{:.1f}', unit = '\\AA')
tab = Table((datecol, pcol, ewcol))
#tab.write(os.path.join(working_dir_path, 'EWtab.tex'), latexdict = ascii.latexdicts['AA'])

#tab.write(os.path.join(working_dir_path, 'EWtab.py'), latexdict = ascii.latexdicts['AA'])

x = w2vsini(wcaII, 393.366 * u.nm).decompose()

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(x, fcaII[0, :], marker = '', drawstyle = 'steps-mid')
ax.set_xlim([-3, +3])
ax.set_xlabel('line shift [v sin(i)]')
ax.set_ylabel('flux')
ax.set_title('Ca II H line in MN Lup (single)')
fig.subplots_adjust(bottom = 0.15)
plt.draw()

#for i in range(len(fcaII)):
 #   print(fcaII.shape[0])
    
#print(fcaII.shape[0]) # rows of fcaII
#print(fcaII.shape[1]) # columns of fcaII
#print(len(fcaII[24]))


#print(type(yshift))
yshift = np.arange((fcaII.shape[0])) * 0.5
#print(yshift)
#print(type(yshift), yshift.dtype)
yshift[:] += 1.5
#print(yshift)
#print(type(yshift), yshift.dtype)
yshift[13:] += 1
#print(yshift)
#print(type(yshift), yshift.dtype)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

#print(fcaII[0])
#print(fcaII[1])
for i in range(25):
    ax.plot(x, fcaII[i,:] + yshift[i], 'k')
    #print(fcaII[i])

#ax.plot(x, np.mean(fcaII, axis = 0))
ax.set_xlim([-2.5, +2.5])
ax.set_xlabel('line shift [$v \\sin i$]')
ax.set_ylabel('flux')
ax.set_title('Ca IIH line in MN Lup (all)')
fig.subplots_adjust(bottom = 0.15)
plt.draw()

'''
 sepearte the grapht into the all and mean
'''
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(x, np.mean(fcaII, axis = 0))
ax.set_xlim([-2.5, +2.5])
ax.set_xlabel('line shift [$v \\sin i$]')
ax.set_ylabel('flux')
ax.set_title('Ca IIH line in MN Lup (mean)')
fig.subplots_adjust(bottom = 0.15)
plt.draw()

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
for i in range(len(fcaII)):
    ax.plot(x, fcaII[i,:] + yshift[i], 'k')
ax.plot(x, np.mean(fcaII, axis = 0))
ax.set_xlim([-2.5, +2.5])
ax.set_xlabel('line shift [$v \\sin i$]')
ax.set_ylabel('flux')
ax.set_title('Ca IIH line in MN Lup (summary)')
fig.subplots_adjust(bottom = 0.15)
plt.draw()



fmean = np.mean(fcaII, axis = 0)
#print(fmean.shape) # (675,)
#print(fcaII.shape) # (25, 675)
#print(fcaII)
#fmean = fmean[np.newaxis, :]
fdiff = fcaII - fmean[np.newaxis, :]
#print(fdiff)
#print('fmean : ', fmean)
#print(fmean.shape) # new shape -> (1, 675) : if fmean is converted by utilizing the 'np.newaxis'
#print(fdiff.shape)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
im = ax.imshow(fdiff, aspect = 'auto', origin = 'lower')

ind1 = delta_p < 1 * u.dimensionless_unscaled
ind2 = delta_p > 1 * u.dimensionless_unscaled
# delta_p : gap of observation time 

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

for ind in [ind1, ind2]:
    im = ax.imshow(fdiff[ind, :], extent = (np.min(x), np.max(x), np.min(delta_p[ind]), np.max(delta_p[ind])), aspect = 'auto', origin = 'lower')
    
ax.set_ylim([np.min(delta_p), np.max(delta_p)])
ax.set_xlim([-1.9, +1.9])
plt.draw()


pplot = delta_p.copy().value
#print(pplot[ind2])
pplot[ind2] -= 1.5
#print(pplot[ind2])
delta_t = np.median(np.diff(delta_p)) / 2.
delta_x = np.median(np.diff(x)) / 2.
fdiff = fdiff / np.max(np.abs(fdiff))

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

for ind in [ind1, ind2]:
    im = ax.imshow(fdiff[ind, :], extent = (np.min(x) - delta_x, np.max(x) + delta_x, np.min(pplot[ind]) - delta_t, np.max(pplot[ind]) + delta_t),
                   aspect = 'auto', origin = 'lower', cmap = plt.cm.Greys_r)
    


ax.set_ylim([np.min(pplot) - delta_t, np.max(pplot) + delta_t])
ax.set_xlim([-1.9, +1.9])
ax.set_xlabel('vel in $v\\sin i$')
ax.xaxis.set_major_locator(plt.MaxNLocator(4))

def pplot(y, pos):
    if y < 0.5:
        yreal = y
    else:
        yreal = y + 1.5
    return yreal

formatter = plt.FuncFormatter(pplot)
ax.yaxis.set_major_formatter(formatter)
ax.set_ylabel('period')
ax.set_title('Final plot in grayscale')
fig.subplots_adjust(left = 0.15, bottom = 0.15, right = 0.99, top = 0.99)
plt.draw()








'''
 check the header of filelist
'''
#f1 = open('C:\Users\Administrator\Desktop\donghun\AA.tex', 'r')
#lines = f1.readlines()
#for line in lines:
#    print(line)
#f1.close()


#for i in range(len(filelist)):
 #   sp1 = fits.open(filelist[i])
  #  head = sp1[0].header
   # sp1.write(os.path.join(working_dir_path, 'test_filelist'), latexdict = ascii.latexdicts['AA'])


'''
sp2 = fits.open(filelist[0])
head = sp2[0].header
print(head)
print(len(sp2)) -> 1
'''