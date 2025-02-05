# emshower.py

import uproot
import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN, HDBSCAN
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.linear_model import RANSACRegressor
from sklearn.linear_model import LinearRegression

# ----------------- Data Processing -----------------

def preprocess_data(df):
    """
    Preprocess the input data for clustering.
    - Extract relevant columns.
    - Sort by event ID, track ID, and Z-coordinate.
    """
    columns_to_keep = ["evtID", "trkID", "segnum", "PID", "X", "Y", "Z", "TX", "TY","E"]
    preprocessed_df = df[columns_to_keep].copy()
    preprocessed_df.sort_values(by=["evtID", "trkID", "Z"], inplace=True)
    return preprocessed_df


# ------ Clustering method: DBSCAN vs HDBSCAN ------

def dbscan_for_pid(df, maxPlate=15, minSample=4, epsilon=15):
    """
    Perform DBSCAN clustering for each PID value and calculate the centroid of each cluster.

    Parameters:
        df (DataFrame): Input DataFrame with columns PID, Z, X, Y, E.

    Returns:
        DataFrame: DataFrame of centroids with columns ['PID', 'Centroid_X', 'Centroid_Y', 'Centroid_Z'].
    """
    centroids = []
    dbscan_labels = pd.Series(index=df.index, dtype=int)  # Initialize labels with the same index as `df`

    # Iterate through each unique PID value
    for pid in df['PID'].unique():
        if pid >= maxPlate:  # Use only the first few plates
            continue

        df_pid = df[df['PID'] == pid].copy()

        # Perform DBSCAN clustering on Z, X, Y coordinates
        # minSample = 2
        if ((pid < minSample) | (len(df_pid)==1)): labels = np.zeros(len(df_pid), dtype=int)
        else: 
            clustering = DBSCAN(eps=epsilon, min_samples=minSample)
            features = df_pid[['X', 'Y']].values
            labels = clustering.fit_predict(features)
        
        dbscan_labels[df_pid.index] = labels  # Collect DBSCAN labels for the entire DataFrame
        # Add DBSCAN labels to the subset DataFrame
        df_pid['labels'] = labels
        # print(pid, df_pid["labels"])
        # print(pid, labels)
        
        # Find the centroids of each cluster (ignoring noise label -1)
        unique_labels = set(labels)
        for label in unique_labels:
            if label > -1:  # Ignore noise points (-1)
                if pid < minSample:
                    cluster_points = df_pid[df_pid['labels'] == label]  # Filter cluster points
                    # sumE = np.sum(cluster_points['E'])
                    # centroid_x = np.sum(cluster_points['E'].values * cluster_points['X'].values) / sumE
                    # centroid_y = np.sum(cluster_points['E'].values * cluster_points['Y'].values) / sumE
                    # centroid_z = np.mean(cluster_points['Z'].values)
                    centroid_x = np.mean(cluster_points['X'].values)
                    centroid_y = np.mean(cluster_points['Y'].values)
                    centroid_z = np.mean(cluster_points['Z'].values)
                else:
                    cluster_points = df_pid[df_pid['labels'] == label]  # Filter cluster points
                    sumE = np.sum(cluster_points['segnum'])
                    XoverE = np.sum(cluster_points['segnum'].values * cluster_points['X'].values)
                    YoverE = np.sum(cluster_points['segnum'].values * cluster_points['Y'].values)
                    centroid_x = XoverE / sumE
                    centroid_y = YoverE / sumE
                    centroid_z = np.mean(cluster_points['Z'].values)

                centroids.append((pid, centroid_x, centroid_y, centroid_z))

    # Convert the centroids list into a DataFrame
    centroid_df = pd.DataFrame(centroids, columns=['PID', 'Centroid_X', 'Centroid_Y', 'Centroid_Z'])

    # Add DBSCAN labels back to the main DataFrame
    df['labels'] = dbscan_labels

    return centroid_df


def hdbscan_for_pid(df, maxPlate=15, minSamples=4):
    """
    Perform HDBSCAN clustering for each PID value and calculate the centroid of each cluster.

    Parameters:
        df (DataFrame): Input DataFrame with columns PID, Z, X, Y, E.
        maxPlate (int): Maximum PID value to consider for clustering.
        minClusterSize (int): Minimum cluster size for HDBSCAN.
        clusterSelectionEpsilon (float): Cluster selection epsilon for HDBSCAN.

    Returns:
        DataFrame: DataFrame of centroids with columns ['PID', 'Centroid_X', 'Centroid_Y', 'Centroid_Z'].
    """

    centroids = []
    hdbscan_labels = pd.Series(index=df.index, dtype=int)  # Initialize labels with the same index as `df`

    # Iterate through each unique PID value
    for pid in df['PID'].unique():
        if pid >= maxPlate:  # Use only the first few plates
            continue

        df_pid = df[df['PID'] == pid].copy()

        # Perform HDBSCAN clustering on Z, X, Y coordinates
        if ((pid < minSamples) | (len(df_pid)==1)): labels = np.zeros(len(df_pid), dtype=int)
        else: 
            clustering = HDBSCAN(min_cluster_size=2, min_samples=2)
        
            features = df_pid[['X', 'Y']].values
            labels = clustering.fit_predict(features)
            # clustering.condensed_tree_.plot(select_clusters=True)

        hdbscan_labels[df_pid.index] = labels  # Collect HDBSCAN labels for the entire DataFrame

        # Add HDBSCAN labels to the subset DataFrame
        df_pid['labels'] = labels
        # print(pid, labels)

        # Find the centroids of each cluster (ignoring noise label -1)
        unique_labels = set(labels)
        for label in unique_labels:
            if label > -1:  # Ignore noise points (-1)
                if pid < minSamples:
                    cluster_points = df_pid[df_pid['labels'] == label]  # Filter cluster points
                    # sumE = np.sum(cluster_points['E'])
                    # centroid_x = np.sum(cluster_points['E'].values * cluster_points['X'].values) / sumE
                    # centroid_y = np.sum(cluster_points['E'].values * cluster_points['Y'].values) / sumE
                    # centroid_z = np.mean(cluster_points['Z'].values)
                    centroid_x = np.mean(cluster_points['X'].values)
                    centroid_y = np.mean(cluster_points['Y'].values)
                    centroid_z = np.mean(cluster_points['Z'].values)
                else:
                    # if label > 0: continue
                    cluster_points = df_pid[df_pid['labels'] == label]  # Filter cluster points
                    sumE = np.sum(cluster_points['segnum'])
                    XoverE = np.sum(cluster_points['segnum'].values * cluster_points['X'].values)
                    YoverE = np.sum(cluster_points['segnum'].values * cluster_points['Y'].values)
                    centroid_x = XoverE / sumE
                    centroid_y = YoverE / sumE
                    centroid_z = np.mean(cluster_points['Z'].values)


                centroids.append((pid, centroid_x, centroid_y, centroid_z))

    # Convert the centroids list into a DataFrame
    centroid_df = pd.DataFrame(centroids, columns=['PID', 'Centroid_X', 'Centroid_Y', 'Centroid_Z'])

    # Add HDBSCAN labels back to the main DataFrame
    df['labels'] = hdbscan_labels

    return centroid_df


# ------ Axis Fitting Method ------

def fit_axis(centroid_df):
    # Fit line in Z-X and Z-Y
    # Z-X fitting
    slope_zx, intercept_zx, r_value_zx, p_value_zx, std_err_zx = stats.linregress(centroid_df['Centroid_Z'], centroid_df['Centroid_X'])
    # Z-Y fitting
    slope_zy, intercept_zy, r_value_zy, p_value_zy, std_err_zy = stats.linregress(centroid_df['Centroid_Z'], centroid_df['Centroid_Y'])

    return (slope_zx, intercept_zx), (slope_zy, intercept_zy)


def fit_robust_axis(centroids, residual_threshold=50):
    """
    Fit a robust axis using centroids and RANSAC.
    
    Parameters:
        centroids (DataFrame): DataFrame of centroids with columns Z, X, Y.
    
    Returns:
        x_model, y_model: Fitted RANSAC models for X-Z and Y-Z.
    """
    z = centroids['Centroid_Z'].values.reshape(-1, 1)
    x = centroids['Centroid_X'].values
    y = centroids['Centroid_Y'].values

    # Fit Z-X relationship using RANSAC
    ransac_x = RANSACRegressor(residual_threshold=residual_threshold)
    try:
        ransac_x.fit(z, x)
        # Check if inliers were found
        if np.sum(ransac_x.inlier_mask_) == 0:
            print("No inliers found for Z-X fitting.")
            slope_zx, intercept_zx = None, None
        else:
            slope_zx = ransac_x.estimator_.coef_[0]  # Extract slope
            intercept_zx = ransac_x.estimator_.intercept_
    except Exception as e:
        print(f"Error during RANSAC Z-X fitting: {e}")
        slope_zx, intercept_zx = None, None

    # Fit Z-Y relationship using RANSAC
    ransac_y = RANSACRegressor(residual_threshold=50)
    try:
        ransac_y.fit(z, y)
        # Check if inliers were found
        if np.sum(ransac_y.inlier_mask_) == 0:
            print("No inliers found for Z-Y fitting.")
            slope_zy, intercept_zy = None, None
        else:
            slope_zy = ransac_y.estimator_.coef_[0]
            intercept_zy = ransac_y.estimator_.intercept_
    except Exception as e:
        print(f"Error during RANSAC Z-Y fitting: {e}")
        slope_zy, intercept_zy = None, None

    return (slope_zx, intercept_zx), (slope_zy, intercept_zy)

# ------ Draw Functions ------

# here df should be the sorted df_filtered which only contains one event
def plot_centroids(df, centroid_df, axis_zx, axis_zy):
    plt.figure(figsize=(12, 6))
    event_id = df['evtID'].values[0]

    # Z-X plot
    plt.subplot(1, 2, 1)
    plt.scatter(df['Z'], df['X'], c='grey', label='Data Points', s=10)
    count = 0
    for pid in df['PID'].unique():
        pid_data = df[df['PID'] == pid]
        clustered = pid_data[pid_data['labels'] > -1]  # Points in a cluster
        non_clustered = pid_data[pid_data['labels'] == -1]  # Noise points

        # Mark clustered points in red and noise points in blue
        
        if(count==0):
            plt.scatter(clustered['Z'], clustered['X'], c='blue', label='cluster', s=10)
            plt.scatter(non_clustered['Z'], non_clustered['X'], c='black', label='noise', s=10)
            count += 1
        else:
            plt.scatter(clustered['Z'], clustered['X'], c='blue',s=10)
            plt.scatter(non_clustered['Z'], non_clustered['X'], c='black',s=10)
        

    plt.scatter(centroid_df['Centroid_Z'], centroid_df['Centroid_X'], c='red', label='Centroids', s=50, marker='x')
    # Plot the fitted axis
    z_values = np.linspace(min(df['Z']), max(df['Z']), 100)
    x_values = axis_zx[0] * z_values + axis_zx[1]
    plt.plot(z_values, x_values, c='green', label='Fitted Axis (Z-X)', linewidth=2)
    plt.plot([], [], ' ', label=f'Event ID: {event_id}')
    plt.xlabel('Z')
    plt.ylabel('X')
    plt.title('Z-X Plot')
    plt.legend()

    # Z-Y plot
    plt.subplot(1, 2, 2)
    plt.scatter(df['Z'], df['Y'], c='grey', label='Data Points', s=10)
    count = 0
    for pid in df['PID'].unique():
        pid_data = df[df['PID'] == pid]
        clustered = pid_data[pid_data['labels'] >= 0]  # Points in a cluster
        non_clustered = pid_data[pid_data['labels'] == -1]  # Noise points

        # Mark clustered points in red and noise points in blue
        if(count==0):
            plt.scatter(clustered['Z'], clustered['Y'], c='blue', label='cluster', s=10)
            plt.scatter(non_clustered['Z'], non_clustered['Y'], c='black', label='noise', s=10)
            count += 1
        else:
            plt.scatter(clustered['Z'], clustered['Y'], c='blue',s=10)
            plt.scatter(non_clustered['Z'], non_clustered['Y'], c='black',s=10)
    plt.scatter(centroid_df['Centroid_Z'], centroid_df['Centroid_Y'], c='red', label='Centroids', s=50, marker='x')
    # Plot the fitted axis
    y_values = axis_zy[0] * z_values + axis_zy[1]
    plt.plot(z_values, y_values, c='green', label='Fitted Axis (Z-Y)', linewidth=2)
    plt.plot([], [], ' ', label=f'Event ID: {event_id}')
    plt.xlabel('Z')
    plt.ylabel('Y')
    plt.title('Z-Y Plot')
    plt.legend()

    plt.tight_layout()
    plt.show()
