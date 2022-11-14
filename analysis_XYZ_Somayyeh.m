% reset the workspace
clear
close all

% load spiral drawing data
d = read_trc("lue-spiral.trc");

% set plotting parameters
TL = [0 5];
nr = 2;
nc = 3;

% plot the left hand marker in x-y-z
marker_name = "L.Finger3.M3";
marker_xyz = d{:,find(names(d) == "L.Finger3.M3") + (0:2)};

t = d{:,"Time"};
t_inds = t>min(TL)&t<max(TL);
t_secs = rem(t(t_inds),1)==0;

% Plot the X, Y, Z data 
figure
subplot(nr,nc,1)
plot (t, marker_xyz, 'LineWidth',2)
legend('X','Y','Z')
title('Raw Data')
ylabel ('mm','FontSize',10)
xlabel ('seconds','FontSize',10)
hold on

% Plot the Y-Z front view
subplot(nr,nc,2)
plot (marker_xyz(:,2), marker_xyz(:,3), 'LineWidth',2)
title('Front View')
ylabel ('Z', 'FontSize', 10)
xlabel ('Y', 'FontSize', 10)
hold on

%%% YOUR CODE HERE
% Filter out large, slow movements with a high-pass butterworth filter at 2
% Hz cutoff and filter out jitter with a low-pass butterworth filter at 20
% Hz cutoff. A 6th order filter is fine.
%%%

% sampling freq fs is the reciprocal of the difference between two points
fs = 1/mean(diff(t));

% cutoff frequencies for the filter
fc_hi = 2;
fc_lo = 20;

% A high-pass butterworth filter at 2 Hz cutoff and a 6th order
[b,a] = butter(6, fc_hi /(fs/2),'high'); 
filtered_marker_xyz_1 = filtfilt(b, a, marker_xyz);

% A low-pass butterworth filter at 20 Hz cutoff and a 6th order
[b,a] = butter(6, fc_lo /(fs/2),'low'); 
filtered_marker_xyz_2 = filtfilt(b, a, filtered_marker_xyz_1);

% A low-pass butterworth filter at 20 Hz cutoff and a 6th order
[b,a] = butter(6, fc_lo /(fs/2),'low'); 
filtered_marker_xyz_3 = filtfilt(b, a, marker_xyz);

% Plot Low Frequency Component
subplot(nr,nc,3)
plot (filtered_marker_xyz_3(:,2), filtered_marker_xyz_3(:,3),'LineWidth',2)
title('Low Frequency Component')
ylabel ('Z', 'FontSize', 10)
xlabel ('Y', 'FontSize', 10)
hold on

% calculate the first PC
[coeff, score, latent, tsquared, explained, mu]=pca(filtered_marker_xyz_2);
pca1_x= score(:,1);
pca1_y= score(:,2);
pca1_z= score(:,3);

subplot(nr,nc,4)
plot (pca1_y,pca1_z, 'LineWidth',1)

hold on
p = polyfit(pca1_y,pca1_z,1); 
f = polyval(p,[-50 50]);
plot([-50 50],f,'r','LineWidth',2) 

title('High Frequency Component and 1st PC')
ylabel ('Z', 'FontSize', 10)
xlabel ('Y', 'FontSize', 10)
hold on
ylim ([-50 50])
xlim ([-50 50])


% calculate projection onto first PC
subplot(nr,nc,4)
proj = filtered_marker_xyz_2*coeff(:,1);


% smooth with a savitsky-golay smoother
proj_smooth = smoothdata(proj,'sgolay');

% count zero crossings
zcd = dsp.ZeroCrossingDetector();
numZeroCross = cast(zcd(proj_smooth(t_inds)),"double");
tremorFrequency = (numZeroCross/2)/max(TL);

% get envelope from 25 sample moving average
env_width = 25;
env = movmax(proj_smooth(t_inds),env_width);

% use the median of the moving maximum as the estimator of the amplitude
amp = median(env);

ttl = round(tremorFrequency,1) + " Hz, " + round(2*amp,1) + " mm amplitude";

% plot
subplot(nr,nc,[5 6])
hold on
plot(t,proj,'k.')
plot(t,proj_smooth,'r')
h1 = refline(0,amp);
h2 = refline(0,-amp);
h1.Color = 0.5*[1 1 1];
h2.Color = 0.5*[1 1 1];
xlim(TL)
title(ttl)
ylabel("mm")
xlabel("seconds")












