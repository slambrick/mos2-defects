% m1=1-0.0099169-0.00032347;
% 
% m2=0-0.0069672-0.00038994;
% m1=1-0.0099169-0.0019826-0.0014433; %as of 27/11/23 nav31
% m2=0-0.0069672-0.00049051-0.00056946; 



% example_point=[-3580 -1080 -750]*10^3; %MOS2 ON HBN
% 
% z_scan2('z_sep', 30e3, 'z_range', [0e3, -6000e3], 'data_path', getappdata(0, 'data_path'),...
% 'example_pos', example_point, 'N_dwell', 250, 'sampling_period', 15,...
% 'make_plots', 1, 'm1', m1, 'm2', m2)



%------------------ multi scans at 200degC monolayer on hbn/sio2
%%multiZ000114 at alpha 130 scan4622

% example_point=[-2696 4.5 -750]*10^3; %mono layer mos2 on hBN
% example_pos_norm_1  = [ -2684 7 -750]*10^3; %hBN
% example_pos_norm_2  = [-2667 -0.5 -750]*10^3; % Silicon substrate
% 
% z_scan_multiple('z_sep', 50e3, 'z_range', [0e3, -6000e3], 'data_path', getappdata(0, 'data_path'),...
% 'example_pos', example_point, 'N_dwell', 250, 'sampling_period', 15,...
% 'make_plots', 1, 'm1', m1, 'm2', m2, 'example_pos_norm_1', example_pos_norm_1, ...
% 'example_pos_norm_2', example_pos_norm_2)

%------------------ pristine mos2 multi scans at 200degC monolayer on hbn/sio2 

m1 = 1-0.0133; %kw570 12/23
m2 = 0-0.0080;
example_point=[-3184 -1225 -750]*10^3; %mono layer mos2 on hBN
example_pos_norm_1  = [ -3250 -1225 -750]*10^3; %sio2
example_pos_norm_2  = [-3186 -1275 -750]*10^3; %bulk mos2


z_scan_multiple('z_sep', 50e3, 'z_range', [0e3, -6000e3], 'data_path', getappdata(0, 'data_path'),...
'example_pos', example_point, 'N_dwell', 250, 'sampling_period', 15,...
'make_plots', 1, 'm1', m1, 'm2', m2, 'example_pos_norm_1', example_pos_norm_1, ...
'example_pos_norm_2', example_pos_norm_2)
