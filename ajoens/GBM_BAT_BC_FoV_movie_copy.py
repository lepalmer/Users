import matplotlib.pylab as plot
import numpy as np
import subprocess


def plot_circle(rad, cx, cy):
    radrad = np.radians(rad)
    r = np.radians(np.arange(0., 365, 5))
    x = np.radians(cx)+radrad*np.cos(r)
    y = radrad*np.sin(r)-np.radians(cy)
    plot.plot(x, y)


def plot_box(rx, ry, cx, cy):
    xarr = np.arange(cx-rx/2., cx+rx/2., 1)
    nx = len(xarr)
    yarr = np.arange(cy-ry/2., cy+ry/2., 1)
    ny = len(yarr)
    x = np.append(np.append(xarr, np.repeat(cx+rx/2., ny)),
                  np.append(xarr[::-1], np.repeat(cx-rx/2., ny)))
    y = np.append(np.append(np.repeat(cy-ry/2., nx), yarr),
                  np.append(np.repeat(cy+ry/2., nx), yarr[::-1]))
    x = np.radians(x)
    y = np.radians(y)
    plot.plot(x, y)


def random_orbit(inclination, x0):
    r = np.radians(np.arange(0., 360, 5))
    x = r-np.pi+x0
    y = np.radians(inclination)*np.sin(r)-0.5*np.radians(inclination)
    return np.degrees(x), np.degrees(y)


gbm_x0 = np.random.rand()*2*np.pi-np.pi
bat_x0 = np.random.rand()*2*np.pi-np.pi
bc_x0 = np.random.rand()*2*np.pi-np.pi
xgbm, ygbm = random_orbit(26, gbm_x0)
xbat, ybat = random_orbit(21, bat_x0)
xbc, ybc = random_orbit(55, bc_x0)

for b in range(len(xbc)):
    plot.clf()
    m = plot.subplot(111, projection='aitoff')
    m.axes.xaxis.set_ticklabels([])
    m.grid(True)
    # GBM
    plot_circle(70, xgbm[b], ygbm[b])
    # BAT
    plot_box(90, 50, xbat[b], ybat[b])
    # BurstCube
    plot_circle(60, xbc[b], ybc[b])
    # print(xgbm[b],ygbm[b])
    # print(xbat[b],ybat[b])
    # print(xbc[b],ybc[b])
    plot.savefig('frame{:03d}.png'.format(b))

subprocess.call(['ffmpeg',
                 '-v', 'warning',
                 '-i', 'frame%03d.png',
                 '-vf', 'palettegen',
                 '-y',
                 'palette.png'])

subprocess.call(['ffmpeg',
                 '-v', 'warning',
                 '-i', 'frame%03d.png',
                 '-i', 'palette.png',
                 '-lavfi', 'paletteuse',
                 '-y',
                 'out.gif'])
