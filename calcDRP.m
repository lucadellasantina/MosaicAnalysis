%% DRPcalc: Calculate Density Recovery Profile
% |Copyright 2017, Luca Della Santina|
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% This software is released under the terms of the GPL v3 software license
%
% *This program calculates the Density Recovery Profile (DRP) of a list of cells. 
% DRP is calculated according to Rodiek method and the software is
% a translation in MATLAB of the original WinDRP 1.6.4 written by Thomas Euler 
% when he was a postdoc in Richard Masland's lab.
%
% Each row of the reference points' matrix is the [X,Y,Z] coordinates of 
% each cell. 
%
% The program performs the following main operations:
%
% # Ask user to choose the cell list 
% # Calculate DRP
% # Plot DRP
%
% *Input:*
%
% * Reference points: one matrix containing coordinates of cells
%  
% *Output:*
%
% DRP plot


[tmpName] = listdlg('PromptString','Select points',...
                           'SelectionMode','single', 'ListString',who);
tmpVars=who; %list available variables
tmpPts=evalin('base',char(tmpVars(tmpName))); %retrieve selected experiment data
clear tmpName tmpVars;

tmpPts2 = tmpPts;
xyPxSize = 0.103; % Number of pixel per micron
%
AOI1 = [500,1500,500,1500]; %TODO: fix this with real AOI values X1, X2, Y1, Y2
AOI2 = [0,2000,0,2000]; %TODO: fix this with real AOI values X1, X2, Y1, Y2

IsCorrectForSoma = true; % When true nCelOverlap counts the cells where soma overlaps
IsCorrectForBorder = false; % When true correct for cells sitting at border of AOI
CalcType = 1;  % 1 cell type(autocorrelation), 2 for crosscorrelaton
IDDS = -1;
IDDS2 = -1;
rAnn = 5;                   % size of each annulus in micron
rAnn = rAnn / xyPxSize;     % express annulus size in pixels
nAnn = 20;                  % number of annuli analyzed
rMax = rAnn*nAnn;           % maxium radius explored
nCelAnn = zeros(1, nAnn);
AreaAOI = -1;
DensAOI = -1;
nCelAOI = -1;
nCelOverlap = 0;
nRndAOIs = 10;

MeanSomaDia = 7; % Average soma diameter of cells in micron
MeanSomaDia = MeanSomaDia / xyPxSize; % express soma diameter in pixels
nCelOverlap = 0;
nCelAOI = 0;

%% Crosscorrelate datasets
% Note: Luca
% For each cell in Dataset1 calculate distance from any cell in Dataset2
% Increment amount of cells in each anulus by the number of cells found
% within that distance (needs to be divided by total number or pairs computed 
% in order to compute the density)

for ir = 1 : size(tmpPts,1)
    % Retrieve ir-th point from dataset1
    xr = tmpPts(ir, 1);
    yr = tmpPts(ir, 2);
    
    % Check if point is in AOI of first image otherwise skip
    if (xr < AOI1(1)) || (xr > AOI1(2)) || (yr < AOI1(3)) || (yr > AOI1(4))
       continue;
    end
        
    nCelAOI = nCelAOI +1;

    % Scan second dataset to measure against point ir-th
    for ic = 1 : size(tmpPts2,1)        
        % Do not process same cell
        if ic == ir, continue; end;
        
        % Retrieve point ic-th from dataset 2
        xc = tmpPts2(ic, 1);
        yc = tmpPts2(ic, 2);

        % Assure point ic-th within AOI of dataset 2
        if (xc < AOI2(1)) || (xc > AOI2(2)) || (yc < AOI2(3)) || (yc > AOI2(4))
            continue;
        end
    
        % Calculate 2D distance between points of interest
        z = sqrt(abs(xr-xc)^2 + abs(yr-yc)^2);
        
        if z < MeanSomaDia
            %Cells are overlapping, note...
            nCelOverlap = nCelOverlap+1;
        end
        
        if IsCorrectForSoma
            % Overlapping cells are skip, anuli start with an offset of
            % MeanSomaDia
            z = z - MeanSomaDia;

            % find which annulus the current cell belongs
            iA = floor(z/rAnn) + 1;
            if (z > 0) && (iA < nAnn)
                nCelAnn(iA) = nCelAnn(iA)+1;
            end
        else
            % No correction at all
            % find which annulus the current cell belongs
            iA = floor(z/rAnn) + 1;
            if (iA > 0) && (iA < nAnn)
                nCelAnn(iA) = nCelAnn(iA) + 1;
            end
        end
    end        
end

nCelCtr = nCelOverlap;
nCelOverlap = nCelOverlap /2;


%% Generate one more histogram with DRP_nCellErr cells more per bin
for ic = 1 : nAnn
   nCelAnnP1(ic) = nCelAnn(ic); %TODO this does not work( + nCellErr); 
end

% Count the cells of data set #2 that are in the AOI defined by set #1
nCelAOI2 = 0;
for ic = 1 : size(tmpPts2,1)
    xc = tmpPts2(ic, 1);
    yc = tmpPts2(ic, 2);
    if (xc > AOI1(1)) && (xc < AOI1(2)) && (yc > AOI1(3)) && (yc < AOI1(4))
        nCelAOI2 = nCelAOI2 + 1;
    end
end

AreaAOI = (AOI1(2) - AOI1(1)) * (AOI1(4) - AOI1(3)) / (xyPxSize^2); % Area evaluated
DensAOI = nCelAOI / AreaAOI;       % Density of cells in AOI of dataset 1
DensAOI2 = nCelAOI2 / AreaAOI;     % Density of cells in AOI of dataset 2
DensAOIm = sqrt(DensAOI*DensAOI2); % Median density

MaxRadi = sqrt(1.1547/DensAOI);     % ??
MaxRadi2 = sqrt(1.1547/DensAOI2);   % ??
MaxRadim = sqrt(1.1547/DensAOIm);   % ??

%% Calculate density profile (density of cells for each annulus)

for i = 1: nAnn
    if IsCorrectForSoma
        % First anulus starts in MeanSomaDia distance from the center of the cell
        AreaAnn(i) = pi * rAnn^2 * (2*(i)-1) + pi*2*MeanSomaDia*rAnn;
    else
        % Fist anulus starts at 0, that means is a circle around the center of the cell
        AreaAnn(i) = pi * rAnn^2 * (2*i-1);
        %lCelAnn(i) = sqrt(nCelAOI*nCelAOI2)*DensAOIm*AreaANN(i);
    end
    DensAnn(i) = nCelAnn(i)/(AreaAnn(i)*sqrt(nCelAOI*nCelAOI2));
    
    %+-one-cell-error bars for the histogram
    DensAnnP1(i) = nCelAnnP1(i)/(AreaAnn(i)*sqrt(nCelAOI*(nCelAOI2))); % +nCellErr*nAnn
end

AreaCtr = 0;
DensCtr = 0;
if MeanSomaDia > 0
    AreaCtr = pi * MeanSomaDia^2;
    DensCtr = nCelCtr/(AreaCtr*sqrt(nCelAOI*nCelAOI2));
end

%% Correct histogram for bounds

if IsCorrectForBorder    
    dxAOI = abs(AOI1(2) - AOI1(1));
    dyAOI = abs(AOI1(4) - AOI1(3));
    c1 = 2/pi *(1/dxAOI + 1/dyAOI);
    c2 = 0.312/dxAOI/dyAOI;
    
    for i=1:nAnn
        if IsCorrectForSoma
            r = (i+0.5)*rAnn + MeanSomaDia;
            eA = 1 - c1*r +c2*r*r;
        else
            r = (i+0.5)*rAnn;
            eA = 1 - c1*r + c2*r*r;
        end
        DensAnn(i) = DensAnn(i)/eA;
        DeansAnnP1(i) = DensAnnP1(i)/eA;
    end
    
    r = 0.5*MeanSomaDia;
    eA = 1 - c1*r + c2*r*r;
    DensCtr = DensCtr/eA;
end

%% Calculate effective radius

Skip1st = (nAnn>4);
if Skip1st
    i=1;
    while i < 3
        Skip1st = Skip1st & (DensAnn(i) < DensAOIm);
        i = i +1;
        % Changed from < DensAOI...., but this should be now correct, since
        % values in DensANN were calculated using the same "mean" density
    end
end

i = 1;
DeadVol = 0;
while (i<nAnn) && ( (DensAnn(i)<DensAOIm) || ( (i==1) && Skip1st))
    DeadVol = DeadVol + (DensAOIm - DensAnn(i)) * (2*i-1);
    i = i + 1;
end

if IsCorrectForSoma
    % DeadVol = DeadVol + (DensAOIm-DensCtr) * (2*i-1);
    % Something like this ?????
end
EffRadi = rAnn * sqrt(DeadVol/DensAOIm);

% Calculate other parameters
PckFact = (EffRadi/MaxRadim)^2;

CtrDens = 1/sqrt(AreaAOI*pi*rAnn*rAnn);
%CtrDens = 1000 * sqrt(DensAOI/(nCelAOI*pi))/rAnn %Alternative formula, the
%former shows the independence of density

Reliabi = DensAOI/CtrDens;
Reliabi2 = DensAOI2/CtrDens;
Reliabim = DensAOIm/CtrDens;

% MinrAnn = Reliabi/(DensAOI*sqrt(pi)*sqrt(AreaAOI));
%% Things to calculate and plot 

% Important factors to caluclate
% Critical density = CrtDens
% Regularity Index = NNind
% Ratio: Data/Random [-] = NNInd/RndNNInd
% Density per annulus = DensAnn[i] for i equal every annulus

bar([1:nAnn]*rAnn * xyPxSize, DensAnn *10^6);
xlabel('Distance from reference cell(um)');
ylabel('Cell density (cells/mm^2)');

clear AOI* Area* Calc* Ctr* EffRadi IDDS* IsCorrect* MaxRadi* NNInd Reliab* Skip nCel*
clear DeadVol MeanSomaDia PckFact RndNNInd Skip1st i* nRndAOIs rMax tmp* xc xr
clear yc yr z* ans;

%% Changelog
%
% _*Version 1.0*             created on 2017-10-03 by Luca Della Santina_
%
%  + Translated WinDRP from Delphi to MATLAB programming language