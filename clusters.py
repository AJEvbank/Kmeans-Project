import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# Create a list of data files
files = os.listdir('.')
data_files = []
for file in files:
    if file.startswith('data_'):
        data_files.append(file)

# Read the data points and cluster assignments
ndata = 0
data = []
clusters = []
for file in data_files:
    with open(file, 'r') as infile:
        count = int(infile.readline().strip())
        ndata += count

        for i in range(count):
            parts = infile.readline().strip().split()
            data.append(float(parts[0]))
            data.append(float(parts[1]))
            clusters.append(int(parts[2]))

data = np.reshape(data, (ndata, 2))

# Read the centroids
centroids = []
with open('centroids.txt', 'r') as infile:
    k = int(infile.readline().strip())

    for i in range(k):
        parts = infile.readline().strip().split()
        centroids.append(float(parts[0]))
        centroids.append(float(parts[1]))

centroids = np.reshape(centroids, (k, 2))

# Plot clusters with different colors
cluster_count = len(set(clusters)) // 2
sns.set_style("whitegrid")
cmap = ListedColormap(sns.color_palette("nipy_spectral", cluster_count)).colors
cmap += ListedColormap(sns.color_palette("terrain", cluster_count)).colors
np.random.shuffle(cmap)
cmap = ListedColormap(cmap)

fig, ax = plt.subplots()
ax.scatter(data[:, 0], data[:, 1], c=clusters, cmap=cmap)
ax.plot(centroids[:, 0], centroids[:, 1], 'k*')
fig.savefig('clusters.png')
plt.show()
