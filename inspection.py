import pylab as P
import numpy as N
from matplotlib.collections import RegularPolyCollection

class MultiPointCursor:
    def __init__(self, axes, useblit=False):
        self.axes = list(axes)
        self.canvas = axes[0].figure.canvas

        self.canvas.mpl_connect('motion_notify_event', self.onmove)
        self.canvas.mpl_connect('draw_event', self.clear)

        self.visible = True
        self.useblit = useblit

        self.loc = [(0,0)]
        self.icons = []
        for ax in self.axes:
            icon = RegularPolyCollection(4,
            sizes=(75,), rotation=0.7853981634,
            facecolors=((0,0,0,0),), edgecolors=('black',),
            offsets=self.loc, transOffset=ax.transData)
            self.icons.append(icon)
            ax.add_artist(icon)
        self.background = None
        self.needclear = False

    def clear(self, event):
        for i in self.icons:
            i.set_visible(False)
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.canvas.figure.bbox)

    def onmove(self, event):
        if not self.visible: return
        if event.inaxes is None: return
        if not self.canvas.widgetlock.available(self): return
        self.needclear = True

        for i in self.icons:
            i.set_offsets([event.xdata, event.ydata])
            i.set_visible(self.visible)

        try:
            self.icons[self.axes.index(event.inaxes)].set_visible(False)
        except ValueError: #event.inaxes not one of our axes
            pass
        
        self._update()

    def _update(self):
        if self.useblit:
            if self.background is not None:
                self.canvas.restore_region(self.background)
            for ax,icon in zip(self.axes, self.icons):
                if self.visible:
                    ax.draw_artist(icon)
            self.canvas.blit(self.canvas.figure.bbox)
        else:
            self.canvas.draw_idle()

        return False

class DataCursor(MultiPointCursor):
    def __init__(self, x, y, data, axes, useblit=True):
        MultiPointCursor.__init__(self, axes, useblit=useblit)
        self.x = x.ravel()
        self.y = y.ravel()
        
        self.data = [d.ravel() for d in data]

        self.xstr = ''
        self.ystr = ''
        
        # we'll let the cursor do the toolbarformatting too.
        for ax in self.axes:
            ax.fmt_xdata = self.fmtx
            ax.fmt_ydata = self.fmty

        self.canvas.mpl_connect('key_press_event', self.key_press_event)

    def onmove(self, event):
        """
        we override event.xdata to force it to snap-to nearest data
        item
        """
        xdata = event.xdata
        ydata = event.ydata

        if xdata and ydata:
            if event.inaxes not in self.axes: return
            curdat = self.data[self.axes.index(event.inaxes)]
            ind = ((self.x - xdata)**2 + (self.y - ydata)**2).argmin()
            event.xdata = self.x[ind]
            event.ydata = self.y[ind]
            valstr = '%.2f' % float(curdat[ind])
            
            otherdata = [d for d in self.data if not d is curdat]
            vals = ['%.2f' % float(dat[ind]) for dat in otherdata]
            valstr = '%s (%s)' % (valstr, ', '.join(vals))
                
            self.xstr = '%1.3f' % event.xdata
            self.ystr = '%1.3f Value=%s' % (event.ydata, valstr)
        MultiPointCursor.onmove(self, event)

    def fmtx(self, x):
        return self.xstr

    def fmty(self, y):
        return self.ystr

    def key_press_event(self, event):
        if event.key == 'c':
            self.visible = not self.visible
            self.needclear = False
            self._update()
            self.canvas.draw()

if __name__ == '__main__':
    from numpy.random import randn
    from matplotlib import cm
    MA = N.ma.MaskedArray

    deg_to_rad = N.pi / 180.

    rng = N.arange(0, 30, .25)
    az = N.arange(0, 361, 1.)
    x = N.sin(az[:,None] * deg_to_rad) * rng[None,]
    y = N.cos(az[:,None] * deg_to_rad) * rng[None,]

    dr = (rng[1:] + rng[0:-1]) / 2.
    da = (az[1:] + az[0:-1]) / 2.

    dx = N.sin(da[:,None] * deg_to_rad) * dr[None,]
    dy = N.cos(da[:,None] * deg_to_rad) * dr[None,]

    data = 10.*randn(*dx.shape)
    mask = data<10
    data = MA(data, mask=mask)
    data2 = MA(data * 1.5 + 5, mask=mask)
    data3 = MA(data - data2, mask=mask)
    data4 = MA(data**2, mask=mask)

    fig = P.figure()
    cmap = cm.jet
    cmap.set_bad('w', 1.0)

    ax = list()
    d_list = [data, data2, data3, data4]
    for panel in xrange(1,5):
        if not ax:
            ax.append(P.subplot(2, 2, panel))
            m = P.pcolormesh(x, y, data, shading='flat', cmap=cmap)
            ax1 = ax[0]
        else:
            ax.append(P.subplot(2, 2, panel, sharex=ax1, sharey=ax1))
            m = P.pcolormesh(x, y, d_list[panel-1], shading='flat', cmap=cmap)
        fig.colorbar(m, ax=ax[-1])
        P.axis('equal')
        P.title('Panel %d' % panel)
        panel += 1

    cursor = DataCursor(MA(dx, mask=mask), MA(dy, mask=mask), d_list, ax)
    P.show()

