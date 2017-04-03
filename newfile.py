"""
copyright (c) 2016 Mahmood Hoseini

this program is free software: you can redistribute it and/or modify
it under the terms of the gnu general public license as published by
the free software foundation, either version 3 of the license, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  see the
gnu general public license for more details.

you should have received a copy of the gnu general public license
along with this program.  if not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
import cPickle
import scipy.io
import scipy.stats
import os

from tep.acquire.database.open_utils import open_database
open_database(str(128437433458131), verbose=False)
from tep.acquire.scripts.read_only_database_shortcut import *
from tep.utils.signal.remove_60_hz import remove_60_hz

display_settings = {
        'axis_label_font_size':20,
        'plot_titles_font_size':24,
        'left_fig_padding':0.075,
        'top_fig_padding':0.075,
        'bottom_fig_padding':0.15,
        'right_fig_padding':0.05,
        'between_fig_x_padding':0.02,
        'between_fig_y_padding':0.02,
        'marker_size':2,
        'axis_linewidth':0.5
        }

plot_titles_font_size = display_settings['plot_titles_font_size']
ticklabel_fontsize = display_settings['axis_label_font_size']
top_plot_titles_y = 0.95
bottom_plot_titles_y = 0.5
left_fig_padding = display_settings['left_fig_padding']
top_fig_padding = display_settings['top_fig_padding']
bottom_fig_padding = display_settings['bottom_fig_padding']
right_fig_padding = display_settings['right_fig_padding']
between_fig_x_padding = display_settings['between_fig_x_padding']
between_fig_y_padding = display_settings['between_fig_y_padding']
axis_linewidth= display_settings['axis_linewidth']

skip_points = 10
subtract_60hz_signal = True
channel = 0
samp_freq = 10000.0
max_num_trials = 18
alpha = 0.6
freq_bands = [[1.0, 2.2],
              [2.3, 4.3],
              [4.8, 9.3],
              [10.0, 24],
              [25.0, 49],
              [52.0, 100],
              [1.0, 20.0],
              [20, 100.0]];

savefig = True
filepath = '/home/davidmorton/turtle-electrophysiology-project/tep_382/oscillations_paper/'
cellsname1 = ['wb070314c4', 'wb071415c1', 'wb021915c1', 'wb022515c1', 'wb052214c2', 'wb071415c2', 
              'wb072215c6', 'wb072215c6', 'wb070314c3', 'wb021715c1', 'wb021915c1', 'wb021915c2', 
              'wb071415c1', 'wb072215c1', 'wb072215c7', 'wb072315c1', 'wb072715c4']
cellsname2 = ['wb070314c5', 'wb071415c2', 'wb021915c2', 'wb022515c2', 'wb052214c3', 'wb071415c3', 
              'wb072215c7', 'wb072215c8', 'wb070314c2', 'wb021715c2', 'wb021915c3', 'wb021915c3', 
              'wb071415c3', 'wb072215c2', 'wb072215c8', 'wb072315c2', 'wb072715c5']


R_lst_alltrials = np.empty((0, len(freq_bands), 2))
R_lst = np.zeros((len(cellsname1), len(freq_bands), 2))
for column, cellname1 in enumerate(cellsname1):
    cellname2 = cellsname2[column]
        
    wavelet_dict = {}
    
    myfile = open(filepath + cellname1 + '.cp', 'rb')
    data1 = cPickle.load(myfile); myfile.close()
    num_trials = len(data1[cellname1]['vr_ids'])
    vr_ids1 = data1[cellname1]['vr_ids']
    channel1 = data1[cellname1]['channel']
    
    myfile = open(filepath + cellname2 + '.cp', 'rb')
    data2 = cPickle.load(myfile); myfile.close()
    num_trials = len(data2[cellname2]['vr_ids'])
    vr_ids2 = data2[cellname2]['vr_ids']
    channel2 = data2[cellname2]['channel']
    
    common_vr_ids = list(set(vr_ids1).intersection(vr_ids2))
    print (cellname1, len(common_vr_ids))
    
    R_ = np.zeros((len(common_vr_ids), len(freq_bands), 2))
    for trial_num, vr_id in enumerate(common_vr_ids):
        VR = VisualRecording.with_id(vr_id)
        onset_time_s = VR.stimulus_onset_ms/1e3
        pre_pulse_s = 2 #2.2
        post_pulse_s = 0.4 #2.8
        duration_s = pre_pulse_s + post_pulse_s
        vtrace1 = VR.voltage_traces_mv[channel1][(onset_time_s-pre_pulse_s)*samp_freq : 
                                                (onset_time_s+post_pulse_s)*samp_freq]
        vtrace2 = VR.voltage_traces_mv[channel2][(onset_time_s-pre_pulse_s)*samp_freq : 
                                                (onset_time_s+post_pulse_s)*samp_freq]
                                               
        if subtract_60hz_signal:
            vtrace1 = remove_60_hz(vtrace1, samp_freq, segment_length=len(vtrace1))
            vtrace2 = remove_60_hz(vtrace2, samp_freq, segment_length=len(vtrace2))
        
        vtrace1 = vtrace1[range(0, len(vtrace1), skip_points)]
        vtrace2 = vtrace2[range(0, len(vtrace2), skip_points)]

        ## Running MatLab Wavelet function
        scipy.io.savemat('/home/davidmorton/Desktop/vmtrace.mat', 
                         {'vm1': vtrace1, 'vm2': vtrace2, 'Fs': samp_freq/skip_points});
        os.system('/usr/local/MATLAB/R2015b/bin/matlab -nojvm -nodisplay <wavelet_vm.m> driver.log')
        mlab_dict = scipy.io.loadmat('/home/davidmorton/Desktop/vm_wavelet.mat')
        wavelet_dict[vr_id] = mlab_dict;
        
        ## Average power & Phase concentration in each freq band
        phase_diff = wavelet_dict[vr_id]['ang2'] - wavelet_dict[vr_id]['ang1']
        R = np.cos(phase_diff)
        P = wavelet_dict[vr_id]['period1']
        for index in range(len(freq_bands)) :
            choice = np.where(np.logical_and(1/P > freq_bands[index][0], 
                                             1/P < freq_bands[index][1]))[1]
            R_[trial_num, index, 0] = R[choice, :pre_pulse_s*samp_freq/skip_points].mean()
            R_[trial_num, index, 1] = R[choice, pre_pulse_s*samp_freq/skip_points:].mean()
    
    R_lst[column, :, :] = np.mean(R_, 0)
    R_lst_alltrials = np.vstack((R_lst_alltrials, R_))
        
scipy.io.savemat('/home/davidmorton/Desktop/Vm_R_list_2.4s.mat', {'R_lst':R_lst, 'R_lst_alltrials':R_lst_alltrials})






