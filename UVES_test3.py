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
#print(fcaII)

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
ax = fig.add.subplot(1, 1, 1)
ax.plot(x, fcaII[0, :], marker = '', drawstyle = 'steps-mid')






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