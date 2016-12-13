% % The purpose of this code is to serve as a tutorial as well as an implementation of a 
% calculation of lagged coherence. Fransen et al. 2016 describe this procedure mathematically, 
% and it is implemented here on simulated LFP data. The point of this tutorial is so explore how 
% transient rhythms can produce low power in individual trials, but create
% spurious high power when averaged across trials. This can arise because Fourier power depends on 
% both amplitude and phase. A high amplitude transient can appear at the same Fourier power peak as a 
% low amplitude rhythmic (persistent) oscillation. When averaged across trials, the power at this 
% frequency can appear continuous (see Jones 2015). One way to observe the typical duration of a rhythm 
% of interest, is to measure the lagged coherence. This measure assesses the ?rhythmicity?. It is a 
% modified calculation of coherence, applied by taking the lagged autospectrum along a trial length to 
% gauge the persistence of rhythms over the trial. It is windowed and tapered in a frequency dependent 
% manner, and is normalized to the amplitude. The result is an index between 0 and 1, indicating the 
% coherence between adjacent epochs along a trial. Choosing longer lags for comparison can inform the 
% reader about the persistence of rhythms of interest within a real data set. A quantitative threshold 
% can be set, guiding choices in data analysis trial length. This procedure will inform the reader about
% the possible transient nature of certain rhythms and allow better interpretation of other trial averaged 
% coherence measures. 

clearvars;
close all;

%% Simulated LFP 
N=7500; %number of time steps
dt=0.001; % sampling interval
T=N*dt;
time_axis=0.001:dt:N*dt; %time vector
trials=50; % number of trials
LFP=zeros(length(time_axis),trials); % LFP will be simulated for each trial for the length of the time vector
LFP_theta_dropin=zeros(length(time_axis),trials); % LFP band that is not persistent throughout the trial- Theta will drop into the LFP signal in bouts

%alpha
alpha=10; % alpha frequency
alpha_phase=rand*pi; % sets the phase if not randomly set on each trial- see for loop

%beta
beta=14;
beta_phase=rand*pi;

%theta
theta=4;
theta_phase=rand*pi;

%intermittent theta
theta_dropin=zeros(1,length(time_axis));
min_bout=100; % specifies the minimum bout length
max_bout=500; % specfifies the maximum bout length
num_bouts=randi([1 10],1); % number of bouts- 
bout_start=randi([1 length(time_axis)-max_bout],num_bouts,1); % randomly set the bout start index along the time vector
bout_duration=randi([min_bout max_bout],1); % varies the bout duration between the min and max

%gamma
gamma=40;
gamma_phase=rand*pi;

% Uncomment phase init to generate random phase offset on each trials
for trial=1:trials

    %alpha_phase=rand*pi;
    %beta_phase=rand*pi;
    %theta_phase=rand*pi;
    %gamma_phase=rand*pi;
   
    alphaLFP=sin(2*pi*alpha*time_axis+alpha_phase); % Generates sinusoid signal at alpha freq and
    betaLFP=sin(2*pi*beta*time_axis+beta_phase);
    thetaLFP=sin(2*pi*theta*time_axis+theta_phase);
    gammaLFP=sin(2*pi*gamma*time_axis+gamma_phase);
    
    LFP(:,trial)=alphaLFP+betaLFP+thetaLFP+gammaLFP; %sums all the sinusoids to make a perfectly rhythmic simulated LFP with no noise.
    
    num_bouts=randi([1 10],1);
    bout_start=randi([0 length(time_axis)-max_bout],num_bouts,1);
    bout_duration=randi([min_bout max_bout],1);
    
        for i=1:num_bouts; %take theta LFP and drops it into the theta dropin LFP signal at the randomized bout start indices.
            theta_dropin(bout_start(i):bout_start(i)+bout_duration)= ...
            thetaLFP(bout_start(i):bout_start(i)+bout_duration);
        end

   LFP_theta_dropin(:,trial)=alphaLFP+betaLFP+theta_dropin+gammaLFP; % Makes the LFP signal with transient theta signal
  
end

%% Fourier transform, fft, Power Spectrum and visualizing simulated data

fj = repmat([(0:N/2-1)*(1/T), (-N/2:-1)*(1/T)],trials, 1)'; 

%  For each freq compute the Fourier transform and save in "X".
X = zeros(size(fj,1),trials);
X_theta_dropin = zeros(size(fj,1),trials);
pow=zeros(size(fj,1),trials);
pow_theta_dropin=zeros(size(fj,1),trials);
tic %normally takes about 200s, run if you have time to kill. >:} or run with fewer trials.
for t=1:trials
  
        for j=1:size(fj,1)
            X(j,t) = sum(LFP(:,t).* (exp(-2*pi*1i*fj(j,t)*time_axis'))); %compute the fourier coefficients manually
            X_theta_dropin(j,t) = sum(LFP_theta_dropin(:,t).* (exp(-2*pi*1i*fj(j,t)*time_axis')));
            %  Compute the power spectrum.
            pow(j,t) = 2*dt^2/T * (X(j,t).*conj(X(j,t)));
            pow_theta_dropin(j,t)= 2*dt^2/T * (X_theta_dropin(j,t).*conj(X_theta_dropin(j,t)));
        end
toc    
end
%  plot some nice figures
figure(1); hold on
plot(fj(:,1), 10*log10(pow(:,1)),'ok',fj(:,1), 10*log10(pow_theta_dropin(:,1)),'*b',...
fj(:,1), 10*log10(mean(pow_theta_dropin(:,:),2)),'+r', 'MarkerSize', 16)
hold off;
legend('Perfectly Rhythmic','One Trial Transient Theta','Average 50 Trials Transient Theta')
xlim([0 50]);  ylim([-50 10]); set(gca,'Fontsize',20);
xlabel('Freq [Hz]'); ylabel('Power');

% you see the same thing when using the fft-
figure(2)
subplot(3,1,1)
periodogram(LFP_theta_dropin(:,1), [], N, 1/dt); xlim([0 60]);  ylim([-10 10]) %peak of the theta frequency is lower
title('LFP with transient theta')
subplot(3,1,2)
periodogram(mean(LFP_theta_dropin(:,:),2), [], N, 1/dt);
title('LFP with transient theta averaged over trials')
xlim([0 60]);  ylim([-10 10])
subplot(3,1,3)
periodogram(mean(LFP(:,:),2), [], N, 1/dt);
title('perfectly rhythmic LFP averaged over trials')
xlim([0 60]);  ylim([-10 10])

%As an illustration...
figure(3); title('Crude look at the power of different frequencies over time in simulated LFP data')
subplot(2,1,1);
spectrogram(LFP(:,1),'yaxis'); ylabel('Normalized Frequency'); ylim([0 0.15]);
xlabel('Time'); 
subplot(2,1,2); title('simulated LFP data with transient theta')
spectrogram(LFP_theta_dropin(:,1),'yaxis'); ylabel('Normalized Frequency'); ylim([0 0.15]);
xlabel('Time');
%% Section 2- Compute Fourier Coefficients with frequency adapted hanning window
% In this example, I will calculate the lagged coherence for a persistent
% and transient theta rhythm embedded in a perfect LFP signal with 3 other
% perfectly rhythmic frequencies (Ref: Fransen et al. 2016)

% cycles of Rhythm of interest
%These values should help guide choices for trade off of freqency resolution vs temporal resolution of lagged coherence
cycles=3; % if too many cycles of a slow frequency are chosen, T may be too small to temporally resolve coherence across lags
freq=theta;% the frequency of interest 
L=cycles*(length(time_axis)/(freq*length(time_axis)*dt));
Nq=1/dt/2;fo=1/(L*dt);
fj = [(0:L/2-1)*(fo), (-L/2:-1)*(fo)];
tiled_windows=repmat(hann(L),floor(length(time_axis)/L),trials);
windowed_LFP=tiled_windows.*LFP;
windowed_LFP_theta_dropin=tiled_windows.*LFP_theta_dropin;

figure(4) %visualize the windowed and tapered LFP 
hold on;
plot(time_axis,windowed_LFP, 'k',time_axis,tiled_windows, 'r');
ylabel('Amplitude', 'Fontsize', 20); xlabel('Time in Seconds', 'Fontsize', 20);
legend('Windowed and Hanning tapered LFP','taper for 3 cycles of theta')

%% Calculating lagged coherence

Ftapered=zeros(L,length(time_axis)/L);
Ftapered_dropin=zeros(L,length(time_axis)/L);

laggedX=zeros(L, length(time_axis)/L-1);
laggedX_dropin=zeros(L, length(time_axis)/L-1);

Fabsq=zeros(L,length(time_axis)/L-1);
Fabsq1=zeros(L,length(time_axis)/L-1);
Fabsq_dropin=zeros(L,length(time_axis)/L-1);
Fabsq1_dropin=zeros(L,length(time_axis)/L-1);
%Fourier Coefficients of tapered epochs
for epoch=1:length(time_axis)/L %here go 12 chunks of 500 since 3 cycles of theta= 500ms
    Ftapered(:,epoch) = fft(windowed_LFP((1+L*(epoch-1)):L*(epoch),1));% calculates fourier coefficients on a window of the LFP
    Ftapered_dropin(:,epoch)=fft(windowed_LFP_theta_dropin((1+L*(epoch-1)):L*(epoch),1)); 
end

for i=1:length(time_axis)/L-1
    %Numerator of eq. 1 from fransen et al.
    laggedX(:,i)=((Ftapered(:,i).*Ftapered(:,i+1))); %lag of one used here- bump up the number of lags here to compare across windows
    laggedX_dropin(:,i)=((Ftapered_dropin(:,i).*Ftapered_dropin(:,i+1)));
    
    %denominator of eq. 1 from fransen et al.
    Fabsq(:,i)=(abs(Ftapered(:,i))).^2;
    Fabsq1(:,i)=(abs(Ftapered(:,i+1))).^2;
    
    Fabsq_dropin(:,i)=(abs(Ftapered_dropin(:,i))).^2;
    Fabsq1_dropin(:,i)=(abs(Ftapered_dropin(:,i+1))).^2;
end
%lagged coherence- normalized to average amplitude 

cohLFP=abs(sum(laggedX(:,:),2)./(sqrt(sum(Fabsq(:,:),2).*sum(Fabsq1(:,:),2))));
cohdropin=abs(sum(laggedX_dropin(:,:),2)./(sqrt(sum(Fabsq_dropin(:,:),2).*sum(Fabsq1_dropin(:,:),2))));
figure(4)
plot(fj,cohLFP, '*-k',fj,cohdropin,'o-r'); xlim([0 60]);
legend('rhythmic LFP', 'transient LFP')
title('Lagged Coherence of persistent and transient simulated LFP data', 'Fontsize', 20)
