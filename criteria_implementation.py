import os
import pandas as pd
import numpy as np

def check_galaxy_masses(mass):
    crit1 = np.ones(mass.size,dtype=bool)
    crit1[(mass < 70) | (mass > 200)] = 0
    return crit1

def check_cluster_masses(data):
    correct_mass = np.zeros(data['row_id'].size,dtype=bool)
    correct_mass[(data['MDPL__FOF__mass'] < 1.8 * 10**15) & (data['MDPL__FOF__mass'] > 0.6 * 10**15)] = 1

    cluster_xs = data['MDPL__FOF__x'][correct_mass]
    cluster_ys = data['MDPL__FOF__y'][correct_mass]
    cluster_zs = data['MDPL__FOF__z'][correct_mass]
    cluster_coordinates = cluster_xs,cluster_ys,cluster_zs
    return cluster_coordinates

def check_distances(cluster_coords,transposed):
    d = np.linalg.norm(transposed-np.asarray(cluster_coords),axis=1)
    d[(d < 8) | (d > 16)] = 0
    return (np.asarray(d).nonzero()[0]) #returns indexes (within crit1) of those that pass crit4 (for given cluster)

def bulk_calc(id):
    d = pd.eval('((df.X - df.X[id])**2 + (df.Y - df.Y[id])**2 + (df.Z - df.Z[id])**2)**(1/2)') #np.linalg.norm is faster above but slower here
    inrange = (d <= 3.125)
    peculiar_x = df.X_Velocity[inrange]
    peculiar_y = df.Y_Velocity[inrange]
    peculiar_z = df.Z_Velocity[inrange]
    vx = peculiar_x.sum() / peculiar_x.size
    vy = peculiar_y.sum() / peculiar_y.size
    vz = peculiar_z.sum() / peculiar_z.size
    v = (vx**2 + vy**2 + vz**2)**(1/2)
    bulk = np.nanmean(v)
    if (622-155) < bulk < (622+155):
        return 1
    else:
        return 0

def contrast(id,mass,vs,rhobar):
    d = pd.eval('((df.X - df.X[id])**2 + (df.Y - df.Y[id])**2 + (df.Z - df.Z[id])**2)**(1/2)')
    mass[d > 3.125] = 0
    rho = mass.sum() / vs
    density_contrast = (rho - rhobar) / rhobar
    print(density_contrast)
    if -0.2 < density_contrast < 3:
        return 1
    else:
        return 0

os.chdir(os.path.dirname(__file__))
sim = input('Please enter path to pickled data (see readme for more info):')
df = pd.read_pickle(sim)

def main():
    masses = df['Mvir']
    crit1 = check_galaxy_masses(masses)

    clusterdata = pd.read_csv('2022-02-08-13-16-34-4734.csv')
    clusters = check_cluster_masses(clusterdata)

    transposed_df_sub = np.transpose([df.X[crit1],df.Y[crit1],df.Z[crit1]])

    crit4 = np.zeros(crit1.nonzero()[0].size,dtype=bool)
    for cluster in np.transpose(clusters): #iterating over the clusters is (slightly) faster than iterating over galaxies
        nearby = check_distances(cluster,transposed_df_sub)
        crit4[nearby] = 1

    crit4_full = np.copy(crit1)
    crit4_full[crit1] = crit4
    crit4_indexes = np.asarray(crit4_full).nonzero()[0]

    crit2 = np.empty(crit4.nonzero()[0].size,dtype=bool)
    for i in range(crit2.size):
        crit2[i] = bulk_calc(crit4_indexes[i])

    print(crit2)
    print(np.asarray(crit2.nonzero()[0]).size)
        
    crit2_full = np.copy(crit4_full)
    crit2_full[crit4_full] = crit2
    crit2_indexes = np.asarray(crit2_full).nonzero()[0]
    np.savetxt('crit2passes_indexes.txt',crit2_indexes)
    np.save('crit2passes.pkl',crit2_full)

    volume_box = 1000**3
    volume_sphere = (4/3) * np.pi * 3.125**3
    rho_bar = masses.sum() / volume_box

    crit3 = np.empty(crit2.nonzero()[0].size,dtype=bool)
    for i in range(crit3.size):
        crit3[i] = contrast(crit2_indexes[i],np.copy(masses.values),volume_sphere,rho_bar)

    crit3_full = np.copy(crit2_full)
    crit3_full[crit2_full] = crit3
    crit3_indexes = np.asarray(crit3_full).nonzero()[0]
    print(crit3_indexes)
    np.save('crit3passes.pkl',crit3_full)
    np.savetxt('full_indexes.txt',crit3_indexes)

if __name__ == '__main__':
    main()