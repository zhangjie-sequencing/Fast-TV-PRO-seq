import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# ==== average_time.py section: Calculate global average pausing time ====
# Read data from chr_read
with open("chr_read") as f:
    chrM = list(map(float, f.readline().strip().split()))  # Convert to float type
    auto = list(map(float, f.readline().split()))

time_points = np.array([0.5, 2, 8, 32])  # Minutes

# Normalization processing
normalized = np.array(auto) / np.array(chrM)
rel_ratio = normalized / normalized[-1]  # Relative ratio

# Define and fit exponential model
def saturation_model(t, k):
    return 1 - np.exp(-k * t)

params, cov = curve_fit(saturation_model, time_points, rel_ratio, p0=[0.1], bounds=(0, np.inf))
k = params[0]
pausing_time_avg = 1 / k  # Global average pausing time

print(f"Decay rate k: {k:.4f} min^-1")
print(f"Global average pausing time: {pausing_time_avg:.2f} minutes")
print(f"Now pausing time estimation could be found in file: pausing_times")

# ==== Pausing_time_fitting.py section: Process each peak and normalize ====
# Use chrM as normalization factor
norm_factors = chrM

# Store data for all peaks and raw pausing times
peak_data = []
pausing_time_raw_list = []

with open('TVPRO_peak') as fin:
    fin.readline()  # Skip header line
    for line in fin:
        cols = line.strip().split()
        chrom, strand, pos = cols[0], cols[1], cols[2]
        reads = list(map(float, cols[4:8]))  # Read reads at each time point
        
        # Normalization processing
        norm_reads = [r / nf for r, nf in zip(reads, norm_factors)]
        
        # Calculate beta estimates
        beta_estimates = []
        for t, nr in zip(time_points, norm_reads):
            if nr >= 1.0:
                nr = 0.999  # Handle saturation value
            if nr <= 0:
                continue
            try:
                beta = -np.log(1 - nr) / t
                beta_estimates.append(beta)
            except:
                continue
        
        # Calculate raw pausing time
        if not beta_estimates:
            pausing_time_raw = 0.0
        else:
            median_beta = np.median(beta_estimates)
            pausing_time_raw = 1 / median_beta if median_beta != 0 else 0.0
        
        peak_data.append((chrom, strand, pos, pausing_time_raw))
        pausing_time_raw_list.append(pausing_time_raw)

# Calculate average of raw pausing times (including zeros)
mean_raw = np.mean(pausing_time_raw_list) if pausing_time_raw_list else 0.0

# Write normalized results
with open('pausing_times', 'w') as fout:
    for data in peak_data:
        chrom, strand, pos, raw = data
        if mean_raw == 0:
            final_time = 0.0
        else:
            normalized = raw / mean_raw  # Normalize to average of 1
            final_time = normalized * pausing_time_avg  # Scale to global average
        fout.write(f"{chrom}\t{strand}\t{pos}\t{final_time:.6f}\n")