clear all;
close all;
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
beatFreq = 77e9;
% Max Range = 200m
maxRange = 200;
% Range Resolution = 1 m
deltaR = 1;
% Max Velocity = 100 m/s
vMax = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 3e8; %m/s

Bsweep = c / (2 * deltaR);
Tchirp = 5.5 * (maxRange *  2 / c);
slope = Bsweep / Tchirp;
%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
targetPos = 77;
targetSpeed = 30;


%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.


%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq

                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
trainCollumns=zeros(1,length(t));


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for a=1:1:size(t,2)         
    
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    r_t(a) = targetPos+(targetSpeed*t(a));
    trainCollumns(a) = 2*r_t(a)/c; % Time delay 
    
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(a) = cos(2 * pi * (fc * t(a) + slope * (t(a)^2)/2));
    Rx(a) = cos(2 * pi * (fc * (t(a) - trainCollumns(a)) + slope * ((t(a)-trainCollumns(a))^2)/2));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(a) = Tx(a) .* Rx(a);
end

%% RANGE MEASUREMENT


 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
reshapeMix = reshape(Mix,[Nr,Nd]);
 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
sigFFT=fft(reshapeMix,Nr)./Nr;
 % *%TODO* :
% Take the absolute value of FFT output
sigFFT = abs(sigFFT);
 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
sigFFT = sigFFT(1:Nr/2);

%plotting the range
figure ('Name','Range from First FFT')

 % *%TODO* :
 % plot FFT output 
 plot(sigFFT);
xlabel('Range (m)');
ylabel('Amplitude (normalized)');
hold on 
[maxVal,marker] = max(sigFFT);
plot(marker,maxVal,'x');
 
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
xlabel('Range');
ylabel('Speed');
zlabel('Amplitude');
title('Radar Output');
%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
trainRows = 9;
trainCollumns = 7;

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation

Gr=5;
Gd=4;

% *%TODO* :
% offset the threshold by SNR value in dB
offSet=1.4;
% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
noiseLevel = zeros(1,1);


% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.


   % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
   % CFAR


RDM = RDM/max(max(RDM)); % Normalizing



% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
 
for a = trainRows+Gr+1 : (Nr/2)-(trainRows+Gr)
    for b = trainCollumns+Gd+1 : (Nd)-(trainCollumns+Gd)
        %Create a vector to store noise_level for each iteration on training cells
        noise_level = zeros(1,1);
        %Step through each of bins and the surroundings of the CUT
        for c = a - (trainRows+Gr) : a + (trainRows+Gr)
            for d = b-(trainCollumns + Gd) : b + (trainCollumns + Gd)
                %Exclude the Guard cells and CUT cells
                if (abs(a - c) > Gr || abs(b - d) > Gd)
                    %Convert db to power
                    noise_level = noise_level + db2pow(RDM(c,d)); %requires SignalToolbox
                end
            end
        end
        
        %Calculate threshould from noise average then add the offset
        threshold = pow2db(noise_level/(2*(trainCollumns+Gd+1)*2*(trainRows+Gr+1)-(Gr*Gd)-1));
        %Add the SNR to the threshold
        threshold = threshold + offSet;
        %Measure the signal in Cell Under Test(CUT) and compare against
        CUT = RDM(a,b);
        
        if (CUT < threshold)
            RDM(a,b) = 0;
        else
            RDM(a,b) = 1;
        end
        
    end
end

RDM(RDM~=0 & RDM~=1) = 0;







% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,RDM);
colorbar;
xlabel('Range');
ylabel('Speed');
zlabel('Amplitude (normalized)');
title('CFAR Output');
 
 