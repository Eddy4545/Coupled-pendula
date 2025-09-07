#Import and Constants

import math
import matplotlib.pyplot as plt

G = 9.81  # Gravitational acceleration (m/s²)
K = 0.1   # Coupling constant (N/m)

delta_t = 0.005   # Time step (s)
t_max = 60        # Total simulation time (s)

# Initial Pendulum Properties
Length_of_pend_1 = Length_of_pend_2 = 0.3   # Initial pendulum length (m)
mass_of_pend_1 = mass_of_pend_2 = 0.05  # Mass of each pendulum (kg)
angle_in_rad_1 = 0.3   # Initial angle for pendulum 1 (rad)
angle_in_rad_2 = 0.0   # Initial angle for pendulum 2 (rad)
anglular_velocity_pend_1 = 0.0 # Initial angular velocity for pendulum 1
anglular_velocity_pend_2 = 0.0 # Initial angular velocity for pendulum 2


# Data Storage
angle_1_values = [] # Angle values for pendulum one
angle_2_values = [] # Angle values for pendulum 2 
peaks = [] # Peak values of angle 1
peak_times = [] # Time when peaks occur
values_of_time = [] # Time point for simulation
lengths = [] # All tested lengths
beat_periods = [] # Beat periods for each tested length
beat_period = 0.0 # Computed beat period


# Function using Verlet's Algorythm (removes data once 1 experiment is conducted)
# Updates angle values without using velocity
def verlet_algotithm():
    global values_of_time, angle_1_values, angle_2_values, beat_period
    values_of_time = []
    angle_1_values = []
    angle_2_values = []    
    beat_period = 0.0
    
    # Building time array from 0 to t_max using intervals of delta_t
    t = 0
    while t < t_max:
        values_of_time.append(t)
        t = t + delta_t
    
    #setting initial angle values using a step from Euler's method to initialize Verlet's algorythm
    angle_1_values = [angle_in_rad_1, angle_in_rad_1 + anglular_velocity_pend_1 * delta_t]
    angle_2_values = [angle_in_rad_2, angle_in_rad_2 + anglular_velocity_pend_2 * delta_t]

    # for each time step we calculate the angular displacement using gravitational acceleration and coupling constant
    for i in range(1, len(values_of_time) - 1):
        acceleration_1 = ((-G / Length_of_pend_1) * angle_1_values[i]) + ((K / mass_of_pend_1) * (angle_2_values[i] - angle_1_values[i]))
        acceleration_2 = ((-G / Length_of_pend_2) * angle_2_values[i]) + ((K / mass_of_pend_2) * (angle_1_values[i] - angle_2_values[i]))
        angle_1_rad_next = (2 * angle_1_values[i]) - angle_1_values[i - 1] + (acceleration_1 * (delta_t ** 2))
        angle_2_rad_next = (2 * angle_2_values[i]) - angle_2_values[i - 1] + (acceleration_2 * (delta_t ** 2))
        angle_1_values.append(angle_1_rad_next)
        angle_2_values.append(angle_2_rad_next)


# Detects max amplitude values from a sequence
# Identifies values in sequence that are closest to the max amplitude to (track repeating peaks)
def function_to_locate_maximum(seq):
    acc = 0.0001
    maximum_valueues_from_sequence = []
    index_from_max_amp_values = []
    maximum_value = angle_in_rad_1
    for a, n in enumerate(seq):
        if n > maximum_value:
            maximum_value = n
        if (maximum_value - acc) < n < (maximum_value + acc):
            maximum_valueues_from_sequence.append(n)
            index_from_max_amp_values.append(a)
    return maximum_valueues_from_sequence, index_from_max_amp_values


# Finds local maxima in the angular position (marks when the pendulum reaches the top of its swing)
def locate_manual_peaks(theta, time):
    peaks = []
    peak_times = []
    for a in range(1, len(theta) - 1):
        if theta[a] > theta[a - 1] and theta[a] > theta[a + 1]:
            peaks.append(theta[a])
            peak_times.append(time[a])
    _, index_from_max_amp_values = function_to_locate_maximum(peaks)
    return peak_times, peaks, index_from_max_amp_values


# Calculate the time between peaks to estimate average beat period
def approximate_time_period(peak_times, threshold=2):
    if len(peak_times) < 2:
        return None

    # Calculate all time differences between consecutive peaks
    periods = [t2 - t1 for t1, t2 in zip(peak_times, peak_times[1:])]
    
    avg = sum(periods) / len(periods)

    # Filter out periods that deviate too much from the average
    filtered = [p for p in periods if abs(p - avg) < threshold]

    return sum(filtered) / len(filtered) if len(filtered) > 2 else None


# Very slightly increases pendulum length
def modify_the_length(step=0.0005):
    global Length_of_pend_1, Length_of_pend_2
    Length_of_pend_1 = Length_of_pend_1 + step
    Length_of_pend_2 = Length_of_pend_2 + step


# Run one simulation with given length
def vary_length():
    global beat_seq, beat_period, peak_times, peaks
    modify_the_length()
    verlet_algotithm()
    peak_times, peaks, beat_seq = locate_manual_peaks(angle_1_values, values_of_time)
    beat_period = approximate_time_period(beat_seq)
    if beat_period is not None:
        lengths.append(Length_of_pend_1)
        beat_periods.append(beat_period)



# Run multiple simulations until it finds optimal length for our target beat period
def period_optimizer(target, acc=0.01, NSTEPS=2000):
    global beat_period
    for i in range(NSTEPS):
        vary_length()
        if beat_period and abs(beat_period - target) < acc:
            print(f'The Length of {Length_of_pend_1:.4f} m gives a beat period of {beat_period:.3f} seconds!')
            break


# Pendulum angles over time with peaks marked
def plot_the_angles_peaks(time, angle_1_rad, angle_2_rad, peak_times , save=True):
    plt.figure(figsize=(15, 7.5))
    plt.plot(time, angle_1_rad, label="Pendulum 1", color='#FF0000')
    plt.plot(time, angle_2_rad, label="Pendulum 2", color='#0000FF')
    values_of_peaks = [angle_1_rad[time.index(t)] for t in peak_times]
    plt.scatter(peak_times, values_of_peaks, color='purple', label="Peaks (Pendulum 1)", zorder=3)
    plt.xlabel("Time (s)")
    plt.ylabel("Angle (rad)")
    plt.title("Angle versus Time with Peaks")
    plt.legend()
    plt.grid(True)
    if save:
        plt.savefig("angle_vs_time_peaks.png")
    plt.show()

# Beat Period VS Pendulum Length
def plot_beat_period_vs_length(save=True):
    plt.figure(figsize=(15, 7.5))
    plt.plot(lengths, beat_periods, label="Length vs Beat Period", color="#FF00FF")
    plt.xlabel("Length (m)")
    plt.ylabel("Beat Period (s)")
    plt.title("Pendulum Length versus Beat Period")
    plt.grid(True)
    plt.legend()
    if save:
        plt.savefig("beat_vs_length.png")
    plt.show()


verlet_algotithm()
period_optimizer(6, acc=0.005, NSTEPS=2000)
plot_beat_period_vs_length()
plot_the_angles_peaks(values_of_time, angle_1_values, angle_2_values, peak_times)

