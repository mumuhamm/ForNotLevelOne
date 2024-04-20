import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# import Axes3D for 3D visualization
from mpl_toolkits.mplot3d import Axes3D
df = pd.read_csv('/Users/am/Documents/06.Tracking/detectors.csv')
df[['volume_id', 'layer_id', 'module_id', 'cx', 'cy', 'cz']].head()
df.head()
x_min, x_max = df['cx'].min(), df['cx'].max()
y_min, y_max = df['cy'].min(), df['cy'].max()
z_min, z_max = df['cz'].min(), df['cz'].max()

print('x: %10.2f %10.2f' % (x_min, x_max))
print('y: %10.2f %10.2f' % (y_min, y_max))
print('z: %10.2f %10.2f' % (z_min, z_max))

df['xyz'] = df[['cx', 'cy', 'cz']].values.tolist()
print(df[['volume_id', 'layer_id', 'module_id', 'xyz']].head())
groupby = df.groupby('volume_id')['xyz'].apply(list).to_frame()
print(groupby)

"""
fig = plt.figure(figsize=(18, 18))

for k in range(groupby.shape[0]):
    ax = fig.add_subplot(3, 3, k+1, projection='3d')
    ax.set_aspect('equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_zlim(z_min, z_max)
    for (idx, row) in groupby.iloc[:k+1].iterrows():
        xyz = np.array(row['xyz'])
        x, y, z = xyz[:,0], xyz[:,1], xyz[:,2]
        ax.plot(x, y, z, linewidth=0.5)
        ax.text(x[0], y[0], z[0], str(idx), None)
        ax.text(x_min, y_min, z_max, r'$\mathbf{FUW}$ $\mathit{Private}$', transform=plt.gca().transAxes, fontsize=12, ha='center', bbox=dict(facecolor='none', edgecolor='none', boxstyle='square'))

plt.tight_layout(pad=0., w_pad=0., h_pad=0.)
plt.savefig('DetectorVisual.pdf')
plt.savefig('DetectorVisual.png', dpi=300)
plt.show()
"""
for k in range(groupby.shape[0]):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect('equal')
    ax.set_xlabel('x', fontsize=24)
    ax.set_ylabel('y', fontsize=24)
    ax.set_zlabel('z', fontsize=24)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_zlim(z_min, z_max)
    for (idx, row) in groupby.iloc[:k+1].iterrows():
        xyz = np.array(row['xyz'])
        x, y, z = xyz[:,0], xyz[:,1], xyz[:,2]
        ax.plot(x, y, z, linewidth=0.5)
        ax.text(x[0], y[0], z[0], str(idx), None)
        ax.text2D(.90, .90, r'$\mathbf{FUW}$ $\mathit{Private}$', transform=ax.transAxes, fontsize=14, ha='right', va='top', bbox=dict(facecolor='none', edgecolor='none', boxstyle='square'))


    plt.tight_layout(pad=0., w_pad=0., h_pad=0.)
    plt.savefig(f'DetectorVisual_{k}.pdf')
    plt.savefig(f'DetectorVisual_{k}.png', dpi=300)
    