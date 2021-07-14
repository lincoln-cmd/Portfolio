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
print(ew.sum())