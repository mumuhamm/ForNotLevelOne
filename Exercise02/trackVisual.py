import os
import matplotlib.pylab as plt
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import trackml
from scipy.stats import pearsonr
from trackml.dataset import load_event
from pandas.plotting import scatter_matrix
import matplotlib.pyplot as plt
import math as math
# Load in convex hull method
from scipy.stats.stats import pearsonr
from scipy.spatial import ConvexHull
#circle
from scipy import optimize
import seaborn as sns






nMinHits=5
cFirstEvent=1010
cEventDataDir='/Users/am/Documents/06.Tracking/train_1'



def getPath(pDataDir,pEventID) :
    return '%s/event%09d' % (pDataDir, pEventID)
hits, cells, particles, truth = load_event(getPath(cEventDataDir,cFirstEvent))
print('Event information: ', particles.head())

#draw track + hits 
def getTrackParameters(pIndex) : 
    dataFrame = pd.DataFrame(particles)
    
# start with those that have 5 hits 
def getTracks(sampleSize) : 
    dataFrame = pd.DataFrame(particles)
    dataFrame = dataFrame[dataFrame['nhits']>=nMinHits]
    # get unique list of particle IDs 
    particle_IDs = np.random.choice(dataFrame.particle_id.unique(),sampleSize)
    print(particle_IDs)
    dataFrame = pd.DataFrame(truth)
    df_truth = dataFrame[dataFrame['particle_id'].isin(particle_IDs)]
    return df_truth

def getHitsFromTracks(df_truth, sampleSize) : 
    dataFrame = pd.DataFrame(hits)
    df_hits = dataFrame[dataFrame['hit_id'].isin(df_truth.hit_id)]
    return df_hits

def getOtherHits(df_truth, sampleSize) : 
    dataFrame = pd.DataFrame(hits)
    df_hits = dataFrame[dataFrame['hit_id'].isin(df_truth.hit_id)== False]
    return  df_hits.sample(n=sampleSize)

#return truths for a given particle 
def getTruth(pTruths, particleID) :
    dataFrame = pd.DataFrame(pTruths)
    df_t = dataFrame[dataFrame['particle_id'] == particleID]
    return df_t


#return hits in a given volume 
def getHitsForVolume(pHits, pVolumeID) : 
    dataFrame = pd.DataFrame(pHits)
    df_v = dataFrame[dataFrame['volume_id'] == pVolumeID]
    #df_v = df_v[df_v['layer_id'] < 6]
    return df_v

#return hits in a given volume 
def getHitsForVolume_perLayer(pHits, pVolumeID, pLayerID) : 
    dataFrame = pd.DataFrame(pHits)
    df_v = dataFrame[dataFrame['volume_id'] == pVolumeID]
    df_v = df_v[df_v['layer_id'] == pLayerID]
    return df_v

# make things look familiar...
#plots hits in (x,y) [cartesian] and (z,r) coordinate system [cylindrical]


def showHitsForVolume(pHits, pVolumeID):
    df_v = getHitsForVolume(pHits, pVolumeID)
    
    # Now estimate r-coordinate (in the x,y plane)
    r = np.sqrt(df_v.x**2 + df_v.y**2)
    phi = np.arctan2(df_v.y, df_v.x)  # Use arctan2 for handling quadrants correctly

    # Creating a new figure with a specified size
    plt.figure(figsize=(6, 6))  # Width, height in inches
    
    # Plotting x vs. y (transverse plane)
    ax1 = plt.subplot(121)
    ax1.scatter(df_v.x, df_v.y, c='blue', edgecolor='none', s=10, alpha=0.7, label='Transverse view')
    ax1.set_xlabel('x [cm]')
    ax1.set_ylabel('y [cm]')
    ax1.set_title('Hit Positions in Transverse Plane')
    ax1.grid(True)  # Adding grid lines for better readability
    ax1.legend()  # Adding a legend to the plot

    # Plotting z vs. r (longitudinal view)
    ax2 = plt.subplot(122)
    ax2.scatter(df_v.z, r, c='green', edgecolor='none', s=10, alpha=0.7, label='Longitudinal view')
    ax2.set_xlabel('z [cm]')
    ax2.set_ylabel('r [cm]')
    ax2.set_title('Hit Positions in Longitudinal Plane')
    ax2.grid(True)
    ax2.legend()

    # Adjusting layout to prevent overlap and ensure clear visibility
    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25, wspace=0.35)

    plt.show()  # Display the plots
    return plt


def showHitsForVolume_perLayer(pHits, pVolumeID, pLayerID) : 
    df_v = getHitsForVolume_perLayer(pHits,pVolumeID,pLayerID)   
    #now estimate r-coordinate (in x,y plane)
    r = (df_v.x**2 + df_v.y**2)**0.5
    phi = np.arctan(df_v.y/df_v.x)
    plt.figure(1)
    plt.subplot(121)
    plt.plot(df_v.x,df_v.y, 'bs')
    plt.xlabel('x [cm]')
    plt.ylabel('y [cm]')

    plt.subplot(122)
    plt.plot(df_v.z,r, 'bs')
    plt.xlabel('z [cm]')
    plt.ylabel('r [cm]')
    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=1.55, hspace=0.25, wspace=0.35)
    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=1.55, hspace=0.25, wspace=0.35)
    
    return plt

def showHitsForParticle(pTruth,particleID) : 
    df_t = getTruth(pTruth,particleID)
    r = (df_t.tx**2 + df_t.ty**2)**0.5
    plt.figure(1)
    plt.subplot(121)
    plt.plot(df_t.tx,df_t.ty, 'bs')
    plt.xlabel('x [cm]')
    plt.ylabel('y [cm]')
    
    plt.subplot(122)
    plt.plot(df_t.tz,r, 'bs')
    plt.xlabel('z [cm]')
    plt.ylabel('r [cm]')
    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=1.55, hspace=0.25, wspace=0.35)    
    return plt

def draw(x,y) : 
    plt.figure(1)
    plt.plot(x,y, 'bs')
    plt.xlabel('x [cm]')
    plt.ylabel('y [cm]')
    
    return plt 

nTrueTracks=1
nFakeHits=5
dh = pd.DataFrame(hits)
dh = dh[np.fabs(dh['z']) < 1]
d_t = getTracks(nTrueTracks)
d_ht = getHitsFromTracks(d_t,nTrueTracks)
d_hf = getOtherHits(d_t,nFakeHits)
r_ht = np.sqrt(d_ht.x**2 + d_ht.y**2)
d_ht['r'] = r_ht
r_hf = np.sqrt(d_hf.x**2 + d_hf.y**2)
d_hf['r'] = r_hf
d = pd.concat([d_ht, d_hf])
plt.plot(dh.x, dh.y, 'o', color='blue')
plt.text(0.13, 1.02, 
         r'$\mathbf{FUW}$ $\mathit{Private}$', 
         transform=plt.gca().transAxes, 
         fontsize=12, ha='center', 
         bbox=dict(facecolor='none', 
                   edgecolor='none', boxstyle='square'))
plt.xlabel('X Label')  # Replace 'X Label' with your actual label
plt.ylabel('Y Label') 
plt.grid(True) 
plt.show()


wdir = os.getcwd()

hits_cols = "hit_id,x,y,z,volume_id,layer_id,module_id,event_name"
particle_cols = "particle_id,vx,vy,vz,px,py,pz,q,nhits,event_name"
truth_cols = "hit_id,particle_id,tx,ty,tz,tpx,tpy,tpz,weight,event_name"
cells_cols = "hit_id,ch0,ch1,value,event_name"

hits_df = pd.DataFrame(columns = hits_cols.split(","))
particle_df = pd.DataFrame(columns=particle_cols.split(","))
truth_df =  pd.DataFrame(columns = truth_cols.split(","))
cells_df = pd.DataFrame(columns= cells_cols.split(','))

def calc_R(xc, yc):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return np.sqrt((x-xc)**2 + (y-yc)**2)

def f_2(c):
    """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
    Ri = calc_R(*c)
    return Ri - Ri.mean()

x = d['x']
y = d['y']
x_m = np.mean(x)
y_m = np.mean(y)

center_estimate = x_m,y_m
center_2, ier = optimize.leastsq(f_2, center_estimate)

xc_2, yc_2 = center_2
Ri_2       = calc_R(*center_2)
R_2        = Ri_2.mean()
residu_2   = sum((Ri_2 - R_2)**2)
print(center_2)
xC = np.linspace((np.min(x)-0.1*R_2), (np.max(x)+0.1*R_2), 100)
yC = np.linspace((np.min(y)-0.1*R_2), (np.max(y)+0.1*R_2), 100)
X, Y = np.meshgrid(xC,yC)
F = (X-xc_2)**2 + (Y-yc_2)**2 - R_2**2
plt.plot(x, y, 'ok')
plt.show()

hits, cells, particles, truth = load_event(getPath(cEventDataDir,cFirstEvent))
print(hits.head())
print(hits.describe())
particles[(particles['q'] != -1) & (particles['q'] != 1)]
print(truth.head())
print(truth.describe())

track = truth[truth['particle_id'] == 4503737066323968]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for hit_id, hit in track.iterrows():
    ax.scatter(hit.tx, hit.ty, hit.tz)
plt.xlabel('X (mm)')
plt.ylabel('Y (mm)')
plt.clabel('Z (mm)')
plt.legend()
plt.show()

def calc_curvature(data_fr):
    x = data_fr.tx
    y = data_fr.ty
    z = data_fr.tz
    ddx  = np.diff(np.diff(x))
    ddy  = np.diff(np.diff(y))
    ddz  = np.diff(np.diff(z))
#     take the mean curvature (not the sum) to avoid bias 
#     since some particles generate more hits and others less
    return np.sqrt(ddx**2 + ddy**2 + ddz**2).mean() 

df  = pd.merge(hits_df,truth_df,how = 'left', on = ['hit_id','event_name'])
df = df[df['particle_id']!= 0] # drop particle 0 
grouped = df.groupby(['event_name','particle_id'])
curvatures = grouped.apply(calc_curvature)


g = sns.jointplot(data=hits, x='x', y='y',  s=1, height=12)
g.ax_joint.cla()
plt.sca(g.ax_joint)

volumes = hits.volume_id.unique()
for volume in volumes:
    v = hits[hits.volume_id == volume]
    plt.scatter(v.x, v.y, s=3, label='volume {}'.format(volume))

plt.xlabel('X (mm)')
plt.ylabel('Y (mm)')
plt.text(0.13, 1.02, 
         r'$\mathbf{FUW}$ $\mathit{Private}$', 
         transform=plt.gca().transAxes, 
         fontsize=12, ha='center', 
         bbox=dict(facecolor='none', 
                   edgecolor='none', boxstyle='square'))
plt.legend()
plt.show()

gz = sns.jointplot(data=hits, x='z', y='y', s=1, height=12)
gz.ax_joint.cla()
plt.sca(gz.ax_joint)


volumesz = hits.volume_id.unique()
for volume in volumesz:
    v = hits[hits.volume_id == volume]
    print(v.z)
    plt.scatter(v.z, v.y, s=3, label='volume {}'.format(volume))

plt.xlabel('Z (mm)')
plt.ylabel('Y (mm)')
plt.text(0.13, 1.02, 
         r'$\mathbf{FUW}$ $\mathit{Private}$', 
         transform=plt.gca().transAxes, 
         fontsize=12, ha='center', 
         bbox=dict(facecolor='none', 
                   edgecolor='none', boxstyle='square'))
plt.legend()
plt.show()


hits_sample = hits.sample(8000)
sns.set(style="whitegrid")
palette = sns.color_palette("viridis", n_colors=len(hits_sample['volume_id'].unique()))  
aspect_ratio = 8/ 4
pairplot = sns.pairplot(data=hits_sample, 
                        hue='volume_id', 
                        height=1.4, aspect=aspect_ratio, 
                        palette=palette, 
                        plot_kws={'s': 20, 'edgecolor': 'w', 'linewidth': 0.5, 'alpha': 0.6})
#pairplot.fig.suptitle(r'$\mathbf{FUW}$ $\mathit{Private}$', fontsize=14)  
pairplot.fig.suptitle(r'$\mathbf{FUW}$ $\mathit{Private}$', fontsize=14, ha='left', x=0.05)

pairplot.fig.subplots_adjust(top=0.95)  # Adjust subplot to give space for the title
pairplot.fig.text(0.5, 1.02, 
                  r'$\mathbf{FUW}$ $\mathit{Private}$',
                    transform=pairplot.fig.transFigure, 
                    fontsize=12, ha='center', va='bottom', 
                    bbox=dict(facecolor='none', edgecolor='none', boxstyle='square'))
plt.savefig('pairplot.pdf')
plt.show()
