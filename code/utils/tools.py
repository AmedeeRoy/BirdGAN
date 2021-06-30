# import python libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import  splrep, splev
import scipy.stats as st
import scipy.signal


R = 6377726
pi = np.pi

def dist_ortho(lon1, lat1, lon2, lat2):
    a = np.sin((lat1 - lat2)/2*pi/180)**2
    b = np.cos(lat1*pi/180)*np.cos(lat2*pi/180)
    c = np.sin((lon1- lon2)/2* pi/180)**2

    dist = R * 2* np.arcsin( np.sqrt(a + b*c))
    return dist

def cap(lon1, lat1, lon2, lat2):
    # to radians
    lat1 = lat1*pi/180
    lat2 = lat2*pi/180
    lon1 = lon1*pi/180
    lon2 = lon2*pi/180

    delta_lon = lon2-lon1

    a = np.cos(lat1) * np.sin(lat2) - np.sin(lat1)*np.cos(lat2)*np.cos(delta_lon)
    b = np.sin(delta_lon) * np.cos(lat2)

    cap = np.arctan2(b , a)
    cap = cap%(2*pi)

    return cap*180/pi

def format_data(data, colony, scale):

    data['x'] = pi * R * (data['lon'] - colony[0])/180
    data['y'] = pi * R * (data['lat'] - colony[1])/180

    data['x_std'] = data['x']/scale
    data['y_std'] = data['y']/scale

    data['dist_colony'] = dist_ortho(colony[0], colony[1], data.lon, data.lat)
    # data['day'] = [t[-5:] for t in data.trip]

    data_new = pd.DataFrame()
    for t in data.trip.unique():
        traj = data[data.trip == t].copy()
        n = len(traj)

        step = dist_ortho( traj.lon.values[0:(n-1)], traj.lat.values[0:(n-1)], traj.lon.values[1:n], traj.lat.values[1:n])
        c = cap( traj.lon.values[0:(n-1)], traj.lat.values[0:(n-1)], traj.lon.values[1:n], traj.lat.values[1:n])
        direction = [d%360 - 360 if d%360 > 180 else d%360 for d in np.diff(c)]

        traj['step_distance'] = np.append(np.nan, step)
        traj['step_direction'] = np.append([np.nan, np.nan], direction)
        data_new = data_new.append(traj, ignore_index=True)

    return data_new

def padding_data(data, padding):

    traj_all = np.zeros((len(data.trip.unique()), 2, padding))
    i = 0
    for tt in data.trip.unique():
        traj = data[data.trip == tt].copy()
        coord = traj.loc[:,('x_std', 'y_std')].to_numpy()
        coord = coord[np.arange(0, coord.shape[0]),:]

        traj_all[i, 0, 1:coord.shape[0]+1] = coord[:,0]
        traj_all[i, 1, 1:coord.shape[0]+1] = coord[:,1]
        # traj_all[i, 2, :coord.shape[0]] = coord[:,2]
        i += 1

    return traj_all

def format_simulation(x, colony, scale):

    data_fake = pd.DataFrame()

    for i in range(x.shape[0]):
        trip = 'simulation_'+ str(i).rjust( 2, '0')
        n = x.shape[2]

        x_std = x[i,0,:]
        y_std = x[i,1,:]
        # dive = x[i,2,:]

        lon = 180 * x_std*scale /pi/R + colony[0]
        lat = 180 * y_std*scale /pi/R + colony[1]

        rows = pd.DataFrame( {'trip':trip, 'lon':lon.tolist(), 'lat':lat.tolist()})
        data_fake = data_fake.append(rows, ignore_index=True)
        data_fake = format_data(data_fake, colony, scale)

    return data_fake

dicolour = { 'blue':   '#1f77b4',  # muted blue
            'orange': '#ff7f0e',  # safety orange
            'green':  '#2ca02c',  # cooked asparagus green
            'red':    '#d62728',  # brick red
            'purple': '#9467bd',  # muted purple
            'brown':  '#8c564b',  # chestnut brown
            'pink':   '#e377c2',  # raspberry     yogurt pink
            'gray':   '#7f7f7f',  # middle gray
          'yellow': '#bcbd22'   # curry yellow-green
          }

def get_trip_duration(data):
    var = []
    for tt in data.trip.unique():
        traj = data[data.trip == tt].copy()
        var.append(len(traj.dist_colony))
    return np.array(var)

def get_trip_dist(data):
    var = []
    for tt in data.trip.unique():
        traj = data[data.trip == tt].copy()
        var.append( np.sum(traj.step_distance) )
    return np.array(var)

# def get_nb_dive(data):
#     var = []
#     for tt in data.trip.unique():
#         traj = data[data.trip == tt].copy()
#         var.append( np.sum(traj.dive > 0.5) )
#     return np.array(var)

def get_trip_sinuosity(data):
    var = []
    for tt in data.trip.unique():
        traj = data[data.trip == tt].copy()
        var.append(2*(np.max(traj.dist_colony) - np.min(traj.dist_colony)) /np.sum(traj.step_distance))
    return np.array(var)


def kde1d(z, bw=2):
    kernel = st.gaussian_kde(z)
    kernel.set_bandwidth(bw_method=bw)
    x = np.arange(min(z), max(z), 1)
    return kernel

def kde2d(data, bw=0.5):
    xx, yy = np.mgrid[-2:2:0.05, -2:2:0.05]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([data.x_std, data.y_std])
    kernel = st.gaussian_kde(values)
    kernel.set_bandwidth(bw_method=bw)
    foraging_area = np.reshape(kernel(positions).T, xx.shape)
    return xx, yy, foraging_area

def get_distribution(data, data_fake):

    plt.figure(figsize = (15,10))

    plt.subplot(2,3,1)
    ## DAILY DISTANCE
    y = get_trip_dist(data)
    y_false = get_trip_dist(data_fake)

    dy = kde1d(y/1e3, bw = 0.5)
    dy_false = kde1d(y_false/1e3, bw = 0.5)

    bins = np.arange(0, 300, 10)
    # plt.hist(y/1e3, alpha = 0.5, edgecolor='k', bins = bins, density = True, color = dicolour['blue'], label = 'True')
    # plt.hist(y_false/1e3, alpha = 0.5, edgecolor='k', bins = bins, density = True, color = dicolour['orange'], label = 'Simulation')
    x = np.arange(0,300, 1)
    plt.hist(y/1e3, color = dicolour['blue'], label = 'True', alpha = 0.5, edgecolor='k', density = True)
    plt.plot(x, dy_false(x), color = dicolour['orange'], label = 'Simulation')
    plt.legend()
    plt.title('Daily Distance (km)')

    score_daily_dist = np.sum((dy_false(x) - dy(x))**2) / np.sum( dy(x)**2)

    plt.subplot(2,3,2)
    ## FORAGING DURATION
    y = get_trip_duration(data)
    y_false = get_trip_duration(data_fake)

    dy = kde1d(y, bw = 0.5)
    dy_false = kde1d(y_false, bw = 0.5)

    bins = np.arange(0, 800, 15)
    # plt.hist(y, alpha = 0.5, edgecolor='k', bins = bins, density = True, color = dicolour['blue'], label = 'True')
    # plt.hist(y_false, alpha = 0.5, edgecolor='k', bins = bins, density = True, color = dicolour['orange'], label = 'Simulation')

    x = np.arange(0,800, 1)
    plt.hist(y, color = dicolour['blue'], label = 'True', alpha = 0.5, edgecolor='k', density = True)
    plt.plot(x, dy_false(x), color = dicolour['orange'], label = 'Simulation')
    plt.title('Foraging Duration (min)')

    score_foraging_duration = np.sum((dy_false(x) - dy(x))**2) / np.sum( dy(x)**2)

    plt.subplot(2,3,3)
    # # NB OF DIVES
    # y = get_nb_dive(data)
    # dy  =  kde1d(y, bw = 0.5)

    # y_false = get_nb_dive(data_fake)
    # dy_false  = kde1d(y_false, bw = 0.5)

    # x = np.arange(0,50, 1)
    # plt.hist(y, color = dicolour['blue'], label = 'True', alpha = 0.5, edgecolor='k', density = True)
    # plt.plot(x, dy_false(x), color = dicolour['orange'], label = 'Simulation')
    # plt.title('Nb of Dives')

    plt.subplot(2,3,4)
    # SINUOSITY
    y = get_trip_sinuosity(data)
    dy  =  kde1d(y, bw = 0.5)

    y_false = get_trip_sinuosity(data_fake)
    dy_false  = kde1d(y_false, bw = 0.5)

    x = np.arange(0,1, 0.05)
    plt.hist(y, color = dicolour['blue'], label = 'True', alpha = 0.5, edgecolor='k', density = True)
    plt.plot(x, dy_false(x), color = dicolour['orange'], label = 'Simulation')
    plt.title('Sinuosity')

    plt.subplot(2,3,5)
    # STEP LENGTH
    y_false = data_fake.step_distance/1e3
    dy_false  = kde1d(y_false.dropna(), bw = 0.5)

    x = np.arange(0,2, 0.01)
    plt.hist(data.step_distance.dropna()/1e3, color = dicolour['blue'], label = 'True', alpha = 0.5, edgecolor='k', density = True)
    plt.plot(x, dy_false(x), color = dicolour['orange'], label = 'Simulation')
    plt.title('Step Length (km)')

    plt.subplot(2,3,6)
    # STEP DIRECTION

    y_false = data_fake.step_direction
    dy_false  = kde1d(y_false.dropna(), bw = 0.5)

    x = np.arange(-180,180, 10)
    plt.hist(data.step_direction, bins = x, color = dicolour['blue'], label = 'True', alpha = 0.5, edgecolor='k', density = True)
    plt.plot(x, dy_false(x), color = dicolour['orange'], label = 'Simulation')
    plt.title('Step Direction (degree)')
    # score_foraging_times = np.sum((dy_false/np.sum(dy_false) - dy/np.sum(dy))**2) / np.sum( (dy/np.sum(dy))**2)

#     return  [score_daily_dist, score_foraging_duration]

def get_periodogram(traj):
    periodogram_lon = []
    periodogram_lat = []

    for j in range(traj.shape[0]):
        lon = traj[j, 0, :]
        lat = traj[j, 1, :]
        x, lon = scipy.signal.periodogram(lon, scaling = 'spectrum', detrend = False)
        x, lat = scipy.signal.periodogram(lat, scaling = 'spectrum', detrend = False)

        periodogram_lon.append(lon)
        periodogram_lat.append(lat)

    return (x, np.array(periodogram_lon), np.array(periodogram_lat))


def get_score(traj, traj_sim, plot=False):
  x, lon, lat = get_periodogram(traj)
  x, lon_GAN, lat_GAN = get_periodogram(traj_sim)

  d_lon = np.log(np.mean(lon, axis=0)/np.mean(lon_GAN, axis=0))**2
  d_lat = np.log(np.mean(lat, axis=0)/np.mean(lat_GAN, axis=0))**2

  if plot:
    plt.figure(figsize = (12,6))

    plt.subplot(1,2,1)
    plt.loglog(x, np.mean(lon, axis=0), linewidth = 3, linestyle = '--')
    plt.fill_between(x, np.mean(lon, axis=0) - 0.5 * np.std(lon, axis=0) ,
                    np.mean(lon, axis=0) + 0.5 * np.std(lon, axis=0),
                    alpha = 0.2)


    plt.loglog(x, np.mean(lon_GAN, axis=0), linewidth = 2)
    plt.title('Longitude')

    plt.subplot(1,2,2)
    plt.loglog(x, np.mean(lat, axis=0), linewidth = 3, linestyle = '--', label = 'True')
    plt.fill_between(x, np.mean(lat, axis=0) - 0.5 * np.std(lat, axis=0) ,
                    np.mean(lat, axis=0) + 0.5 * np.std(lat, axis=0),
                    alpha = 0.2)

    plt.loglog(x, np.mean(lat_GAN, axis=0), linewidth = 2, label = 'GAN')

    # plt.xlim([0, 0.5])
    # plt.ylim([1e-6, 1e1])
    plt.legend()
    plt.title('Latitude')

    plt.show()

  return np.mean(d_lon + d_lat)/2
