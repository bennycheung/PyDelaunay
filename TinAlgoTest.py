# TinAlgoTest.py
#
# Testing for TinAlgo module, after randomly insert a set of points
# using matplotlib to plot the delaunay triangles and voronoi polygons.
#
import random 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from TinAlgo import *

# Draw the delaunay triangles
def draw_delaunay(builder, sites):
    for site in sites:
        vpList = builder.getSiteDelaunay(site)
        for vp in vpList:
            plt.plot([site.x,vp.x], [site.y,vp.y], 'r--')

# Draw the voronoi polygons
def draw_voronoi(builder, sites):
    for site in sites:
        vpList = builder.getSiteVoronoi(site)
        vpList.append(vpList[0])
        plt.plot([p.x for p in vpList], [p.y for p in vpList], 'b', lw=2)

# initialize the TIN builder with the test range of X=[0,100] and Y=[0,100]
builder = TinBuilder(0.0, 100.0, 0.0, 100.0)

# randomly generate a set of points
sites = []
num = 100
i = 0
while i < num:
    site = Site(random.randint(0,100), random.randint(0,100))
    print site.x, site.y
    site.sitenum = i
    sites.append(site)
    edge = builder.insertSite(site)
    if (edge != None):
        i += 1
    
# output the TIN results
draw_delaunay(builder, sites)
draw_voronoi(builder, sites)

plt.plot([p.x for p in sites], [p.y for p in sites], 'ro')
plt.axis([0,100,0,100])
plt.show()

