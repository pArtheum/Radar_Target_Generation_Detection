clear;
close all;
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Radar Carrier Frequency
main_fc = 77e9;
% Radar Max RANGE
max_range = 200;
% Range Resolution
range_res = 1;
% Maximum Velocity
max_vel = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Speed of light
cel = 3e8;

%% User Defined Range and Velocity of target
% define the target's initial position and velocity. Note : Velocity
% remains constant
tgt_Pos = 50;
tgt_vel = -30;


%% FMCW Waveform Generation

%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.
B = cel / (2 * range_res);
Tchirp = 5.5 * ((2 * max_range) / cel);
slope = B / Tchirp;
                                                         
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples

%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %t_cell_rangeansmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td =zeros(1,length(t));


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 
for i=1:length(t)         
    %For each time stamp update the Range of the Target for constant velocity. 
    timestamp = (tgt_Pos + (t(i) * tgt_vel)) / cel;

    %For each time sample we need update the t_cell_rangeansmitted and
    %received signal. 
    Tx(i) = cos(2 * pi * ((main_fc * t(i)) + (slope * t(i)^2) / 2));
    tmp = (t(i) - timestamp);
    Rx(i) = cos(2 * pi * ((main_fc * tmp) + (slope * tmp^2) / 2));
    
    %Now by mixing the t_cell_rangeansmit and Receive generate the beat signal
    %This is done by element wise mat_cell_rangeix multiplication of t_cell_rangeansmit and
    %Receiver Signal
    Mix(i) = Rx(i) * Tx(i);
end

%% RANGE MEASUREMENT

%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
beat = reshape(Mix, [Nr, Nd]);

%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
beat_fft = fft(beat) / Nr;

% Take the absolute value of FFT output
abs_beat_fft = abs(beat_fft);

% Output of FFT is double sided signal, but we are interested in only one side of the spect_cell_rangeum.
% Hence we throw out half of the samples.
side_fft = abs_beat_fft(1 : Nr / 2);

%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)

% plot FFT output 
f = Nr / length(side_fft) ;
plot(f * (0 : (Nr / 2 - 1)), side_fft);
axis ([0 200 0 1]);



%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

%Select the number of t_cell_rangeaining Cells in both the dimensions.
t_cell_range = 8;    % For Range Axis
t_cell_doppler = 2;  % For Doppler Axis

%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
g_cell_range = 2;   % For Range Axis
g_cell_doppler = 1;  % For Doppler Axis

% offset the threshold by SNR value in dB
% After tests, iffset of 10 keep the nois and 12 keep a secondary match 
offset = 14;

%Create a vector to store noise_level for each iteration on t_cell_rangeaining cells
noise_level = zeros(1,1);

tot_guard_cells = (2 * g_cell_range + 1) * (2 * g_cell_doppler + 1) -1;
tot_training_cells = (2 * t_cell_range + 2 * g_cell_range + 1) * (2 * t_cell_doppler + 2 * g_cell_doppler + 1) - (tot_guard_cells + 1);

%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for t_cell_rangeaining and Guard Cells.
%For every iteration sum the signal level within all the t_cell_rangeaining
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the t_cell_rangeaining
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.

CFAR = zeros(size(RDM));
% Use RDM[x,y] as the mat_cell_rangeix from the output of 2D FFT for implementing
% CFAR

% Set window offsets as NESW from cent_cell_rangeal pixel
offset_h = t_cell_range + g_cell_range;
offset_w = t_cell_doppler + g_cell_doppler;

% creating position vector of pixels to process by eliminating window
% offset
range_vector = offset_h + 1 : Nr/2 - offset_h;
doppler_vector = offset_w + 1 : Nd - offset_w;

for range_idx = range_vector
    for doppler_idx = offset_w + 1 : Nd - offset_w
        % Get sliding window around evaluated pixel
        window = RDM(range_idx - offset_h : range_idx + offset_h, ...
                       doppler_idx - offset_w : doppler_idx + offset_w);

        % As the window contains guard cells and evaluated cell, we set
        % them to 0 to not influence threshold evaluation
        window(range_idx - g_cell_range : range_idx + g_cell_range, ...
                 doppler_idx - g_cell_doppler : doppler_idx + g_cell_doppler) = 0;

        % For every iteration sum the signal level within all the training 
        % cells. To sum convert the value from logarithmic to linear using
        % db2pow function.
        window = sum(db2pow(window));

        %nAverage the summed values for all of the training cells used. 
        % After averaging convert it back to logarithmic using pow2db.
        % Revert average power to decibels
        window = pow2db(window / tot_training_cells);

        % Add the offset to determine the SNR threshold
        threshold = window + offset;

        % Apply the threshold to the CUT
        if RDM(range_idx, doppler_idx) > threshold
            CFAR(range_idx, doppler_idx) = 1;
        end
    end
end

%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,CFAR);
colorbar;


 
 