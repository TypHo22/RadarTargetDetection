# RadarTargetDetection
Example code for radar detection with matlab

<b>Using the given system requirements, design
a FMCW waveform. Define its Bandwidth (B), chirp time (Tchirp) and slope of the chirp.</b>
![image](https://user-images.githubusercontent.com/42981587/136083523-32aa7553-848c-4c45-95ee-adc8b489f4f6.png)




<b>Simulate Target movement and calculate the beat or mixed signal for every timestamp.</b>
![image](https://user-images.githubusercontent.com/42981587/136083774-7c5281d0-1da9-4613-ae73-b91904e23b51.png)

FFT for target at 77 meters
![image](https://user-images.githubusercontent.com/42981587/136084048-9de3cdd0-bd14-451a-a72c-4bd57b173ea1.png)

<b>Implement the Range FFT on the Beat or Mixed Signal and plot the result.</b>
Plot show maximum of amplitude normalized over speed and range. Noisy data around peak gets removed with CFAR algorithm. 
![image](https://user-images.githubusercontent.com/42981587/136084500-30bb69fa-3c0b-48a0-adc8-63345d600ee3.png)



<b>Implement the 2D CFAR process on the output of 2D FFT operation, i.e the Range Doppler Map.</b>

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


```matlab
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
```
Filtered Signal:
![image](https://user-images.githubusercontent.com/42981587/136084760-958aa9f7-48e1-4091-80f8-2e7fa1417564.png)

