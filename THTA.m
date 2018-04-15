% THTA is an Electromagnetic Transient Program
%% About the program
% FEDERAL UNIVERSITY OF ITAJUBA - UNIFEI
% INSTITUTE OF ELECTRICAL SYSTEMS AND ENERGY - ISEE
%
% THE UNIVERSITY OF BRITISH COLUMBIA - UBC
% DEPARTMENT OF ELECTRICAL AND COMPUTER ENGINEERING - ECE
%
% Prof. Ph.D. Benedito Donizeti Bonatto - UNIFEI
% Prof. Ph.D. José R. Martí - UBC
% Prof. Dr-Ing. Hermann W. Dommel - UBC
%
% Aluno Luis Fernando Ribeiro Ferreira - UNIFEI
% Aluno Alexandre Fonseca Rocha Barreto - UNIFEI
% Aluno Natanael de Souza Figueiredo - UNIFEI
%
% Copyright (c) 2013 - 2015. All rigths reserved.
%% This version includes:
% - Basics RLC components
% - AC and DC voltage and current sources
% - Triangular current and voltage source
% - Switches controlled by time
% - THTA scheme to elliminate numerical oscillations
%% The input file format
% To use this program you must following structure:
% -Every line should have 10 cloumns
% - The first line MUST BE the Time Card
%     'T' (Number of Nodes in the circuit) (Number of Voltages Sources) (Step Size) (Maximum time for simulation) 0 0 0 0 0
%
%     For instance, the line
%     T   3   1   100E-6  50E-3   0   0   0   0   0
%
%          3 = Number of Nodes in the circuit
%          1 = Number of Voltages Sources
%     100E-6 = Step Size = 100 micro seconds
%      50E-3 = Maximum time for simulation = 50 mili seconds = 3 * 60Hz cycles
%      5 * 0 to fill 10 columns
%
% - The following lines specify the elements in the circuit.
%     Resistor
%     (Element Type) (Node K) (Node M) (Resistance)   0 0 0 0 0 (Plt)
%           R            1        0      20           0 0 0 0 0   3
%
%           R  = Resistor Type
%           1  = Node K
%           0  = Node M
%           20 = Resistance
%           0 0 0 0 0   - To fill 10 columns
%           3  = Plot type 3 - Voltage and current
%
%     Inductor
%     (Element Type) (Node K) (Node M) (Inductance in mH)   0 0 0 0 0 (Plt)
%           L            2        1          100            0 0 0 0 0   5
%
%           L   = Inductor Type
%           2   = Node K
%           1   = Node M
%           100 = Inductance in mH = 100 mH
%           0 0 0 0 0   - To fill 10 columns
%           5  = Plot type 5 - Voltage, current, power and energy
%
%     Capacitor
%     (Element Type) (Node K) (Node M) (Capacitance in uF)    0 0 0 0 0 (Plt)
%           C            3        2          500              0 0 0 0 0   5
%
%           C   = Capacitor Type
%           3   = Node K
%           2   = Node M
%           500 = Capacitance in uF = 500 uF
%           0 0 0 0 0   - To fill 10 columns
%           5  = Plot type 5 - Voltage, current, power and energy
%
%     Continuous voltage source
%     (Element Type) (Node K) (Node M)  (Voltage in V)    0 0 0 0 0 (Plt)
%           EDC          3        2          12           0 0 0 0 0   5
%
%           EDC = Continuous voltage source Type
%           3   = Node K (from)
%           2   = Node M (to)
%           12 = Voltage in V = 12 V
%           0 0 0 0 0   - To fill 10 columns
%           5  = Plot type 5 - Voltage, current, power and energy
%
%     Sinusoidal voltage source
%     (Element Type) (Node K) (Node M) (Peak voltage in V) (Phase angle in deg) (Frequency in Hz) 0 0 0 (Plt)
%           EAC          3        2          191                   90                   60        0 0 0   5
%
%           EAC = Sinusoidal voltage source Type
%           3   = Node K (from)
%           2   = Node M (to)
%           191 = Peak voltage in V = 191 V
%           90  = Initial phase angle in degrees = 90º
%           60  = Frequency in Hz = 60 Hz
%           0 0 0   - To fill 10 columns
%           5  = Plot type 5 - Voltage, current, power and energy
%
%     Triangular voltage source
%     (Element Type) (Node K) (Node M) (Peak voltage in V) (Time of peak voltage) (Time of half voltage) 0 0 0 (Plt)
%           ETRI         3        2          191                   3E-3                   30E-3          0 0 0   5
%
%           ETRI = Triangular voltage source Type
%           3   = Node K (from)
%           2   = Node M (to)
%           191 = Peak voltage in V = 191 V
%           3E-3 = Time of peak voltage = 3 ms
%           30E-3 = Time of half of peak voltage on droping stage = 30 ms
%           0 0 0   - To fill 10 columns
%           5  = Plot type 5 - Voltage, current, power and energy
%
%     Transmission Line
%     (Element Type) (Node K) (Node M) (Characteristic Impedance) (Propagation Time)  0 0 0 0 (Plt)
%           TL           4        3               400                    1E-3         0 0 0 0   5
%
%           TL   = Transmission Line Type
%           4    = Node K
%           3    = Node M
%           400  = Characteristic Impedance Zc in Ohms = 400 Ohms
%           1E-3 = Propagation Time in seconds = 1 mili second
%           0 0 0 0   - To fill 10 columns
%           5  = Plot type 5 - Voltage, current, power and energy
%
%      Switch
%     (Element Type) (Node K) (Node M) (Zero) (TClose) (TOpen) 0 0 0 (Plt)
%           S           5        4        0     1E-3    20E-3  0 0 0   5
%
%
% - The required node voltages may be put in the last line. Use one column
%   for each node. Should you need to plot more than 9 plots, create other
%   'NV' line, as much as needed. Remember to complete the line with 10
%   columns. If no 'NV' is used, the program will ask the nodal voltages to
%   plot at the end of the simulation.
%
% NV 1 2 3 4 0 0 0 0 0
%
%       NV = Nodal Voltage
%       The nodes 1, 2, 3 and 4 will be plotted. The following 0 0 0 0 0
%       are used to complete the 10 columns.
%% Plot the Branches signals of voltage, current, power and energy
% Column 10 in the input data file determines the output plotting:
%  1 -> Current
%  2 -> Voltage
%  3 -> Currente and Voltage
%  4 -> Power and Energy
%  5 -> Current, Voltage, Power and Energy

%% Clear variables
clc;
clear all;
close all;
format short eng;

%% Version and date
% For each new version released, please modify the VersionNumber and
% VersionDate as needed.
VersionNumber = '0.45';
VersionDate   = '10/Sep/2014 17:52';

%% Read the input data
fprintf('UNIVERSIDADE FEDERAL DE ITAJUBÁ \n');
fprintf('THE UNIVERSITY OF BRITISH COLUMBIA \n');
fprintf('Electromagnetic Transients Program \n');
fprintf(['Version ' VersionNumber ', ' VersionDate ' - For academic use \n']);
fprintf('-------------------------------------------------------\n');
FILENAME = input('Enter the input file name *.txt : ', 's');
fprintf('-------------------------------------------------------\r');
fprintf('-------------------- Processing ... -------------------\r');
FID  = fopen(FILENAME);
data = textscan(FID, '%s %d %d %f %f %f %f %d %d %d');
fclose(FID);

% data contains the input file organized as a table. Now, each column
% should be put in a separated column vector.
type = data{1}; from = data{2}; to   = data{3}; val4 = data{4}; val5 = data{5};
val6 = data{6}; val7 = data{7}; val8 = data{8}; val9 = data{9}; plt  = data{10};

% bmax is the number of lines in the input file.
bmax = length(type);


%% Calculating parameters for vector and matrix dimensioning
b = 1;
if strcmp(type(b), 'T')
    % Reading time card data
    % N = Number of nodes
    % M = Number of voltage sources
    % dt = step size
    % tmax = maximum time of simulation
    N = from(b);
    M = to(b);
    dt = val4(b);
    tmax = val5(b);
    Damp = plt(b);
    if Damp == 0
        fprintf('-------------------- THTA Activated --------------------\r');
    end
else
    error('Input file error: the first line is not the Time Card.');
end
PlotNodes = zeros(1, 1);
v = 1;
for b = bmax :-1: 2 % Decreasing loop from bmax until 2
    if strcmp(type(b), 'NV')
        bmax = bmax - 1;
        for c = 2 : 10
            if data{c}(b) > 0
                % Store this node for further plotting
                PlotNodes(v) = data{c}(b);
                v = v + 1;
            else
                break; % No more required plots.
            end
        end
    else
        break; % No more NV lines.
    end
end

% nh = number of energy storage elements
% nTL = number of transmission lines
% nSW = number of switches
% ndcvs = number of DC voltage sources
% nacvs = number of AC voltage sources
% natvs = number of Triangular voltage sources
% ndccs = number of DC current sources
% naccs = number of AC current sources
% natcs = number of Triangular current sources

% number of plots required
nh = 0;
nTL = 0;
nSW = 0;
ndcvs = 0;
nacvs = 0;
ntrvs = 0;
ndccs = 0;
naccs = 0;
ntrcs = 0;

nplt = 0;

for b = 2:bmax
    % Set tstart to zero. If we find any source that was active at t<0,
    % set the program to compute the Steady State Solution.
    tstart = 0;
    if strcmp(type(b), 'L') || strcmp(type(b), 'C')
        % finding the number of lumped elements which need history terms
        nh = nh + 1;
        val8(b) = nh;
    elseif strcmp(type(b), 'TL')
        % finding the number of distributed elements which need history terms
        nTL = nTL + 1;
        val6(b) = floor(val5(b)/dt);
        val7(b) = (val6(b)+1)*dt-val5(b);
        if val6(b) == 0
            error('Delta t must be less than any time of propagation on transmission lines')
        end
        val8(b) = nTL;
    elseif strcmp(type(b), 'S')
        nSW = nSW + 1;
        val8(b) = nSW;
        
        if plt(b) == 0
            % Switches MUST have their instant current calculated.
            nplt = nplt + 1;
            val9(b) = nplt;
        end
        
    elseif strcmp(type(b), 'EDC')
        % finding the number of voltage sources
        ndcvs = ndcvs + 1;
    elseif strcmp(type(b), 'EAC')
        nacvs = nacvs + 1;
    elseif strcmp(type(b), 'ETRI')
        ntrvs = ntrvs + 1;
    elseif strcmp(type(b), 'IDC')
        % finding the number of current sources
        ndccs = ndccs + 1;
    elseif strcmp(type(b), 'IAC')
        naccs = nacvs + 1;
    elseif strcmp(type(b), 'ITRI')
        ntrcs = ntrcs + 1;
    end
    if plt(b) > 0
        % finding the number of elements which have required output plotting
        nplt = nplt + 1;
        val9(b) = nplt;
    end
end
nvs = ndcvs + nacvs + ntrvs;
ncs = ndccs + naccs + ntrcs;

%% Calculate the number of points for the maximum simulation time
npoints = fix(tmax/dt) + 1;

%% Pre-Alocate matrices and vectors
% Recall that
% N = Number of nodes
% M = Number of voltage sources
% Therefore, D is the number of unknown voltage nodes.
D = N-M;

% gkm = conductance for each element in the input file
gkm = zeros(bmax,1);
G = zeros(N,N);
I = zeros(N,1);
IA = zeros(D,1);
V = zeros(N,npoints);
VA = zeros(D,1);
VB = zeros(M,1);

t = zeros(npoints, 1);

% Ih is the history current source for the energy-storage elements.
Ih = zeros(nh,1); Ihnew = Ih; Ihold = Ih;
% Ihkm and Ihmk are the history sources for the transmission lines
Ihkm = zeros(nTL,npoints+1);
Ihmk = zeros(nTL,npoints+1);
IhkmTAU = zeros(1,nTL);
IhmkTAU = zeros(1,nTL);
Ifc = zeros(1,bmax);
% Store the voltage, current, power and energy for each required plot
vkm = zeros(nplt,npoints);
ikm = zeros(nplt,npoints);
imk = zeros(nplt,npoints);
pkm = zeros(nplt,npoints);
pmk = zeros(nplt,npoints);
ekm = zeros(nplt,npoints);
emk = zeros(nplt,npoints);
%Switch control vectors
node_pointer = zeros(1,N);
sw_status = zeros(1,nSW);

%% Initial conditions for switches
% This program represents a switch as a resistor, either a very large or a
% very small one.
% Typical Values:
RswOpen = 1e9;
RswClosed = 1e-6;

SwOper = zeros(nSW, 2);

for b = 2:bmax
    if strcmp(type(b), 'S')
        
        % val8 hold the index to link the switch with the information in
        % the 'SwOper' variable. Column 1 and Columns 2 hold the
        % information if the switch already opened or closed, respectively.
        val5(b) = dt*round(val5(b)/dt);
        val6(b) = dt*round(val6(b)/dt);
        tclose  = val5(b);
        topen   = val6(b);
        idx     = val8(b);
        
        if ( tclose > 0 && topen > 0 )
            % The switch will act after the start of the simulation
            if tclose < topen       % If the switch closes before open
                val4(b) = RswOpen;  % it must be open.
                sw_status(idx) = 0;
            elseif topen < tclose   % If the switch opens before close
                val4(b) = RswClosed;% it must be closed.
                sw_status(idx) = 1;
            end
        end
        
        % Particular cases
        % Initial condition (at t=0) with switch l closed
        if ( tclose < 0 || topen == 0)
            val4(b) = RswClosed;
            sw_status(idx) = 1;
            SwOper(idx, 1) = 1;
        end
        
        % Initial condition (at t=0) with switch l opened
        if ( tclose == 0 || topen < 0)
            val4(b) = RswOpen;
            sw_status(idx) = 0;
            SwOper(idx, 1) = 1;
        end
        
    end
end


%% Build and partition the G matrix
% Calculate branch conductances for each input data row
for b = 2:bmax
    if strcmp(type(b), 'R') || strcmp(type(b), 'S')
        R = val4(b);
        gkm(b) = 1/R;
    elseif strcmp(type(b), 'L')
        L = val4(b)*1e-3; %mH
        gkm(b) = dt/(2*L);
    elseif strcmp(type(b), 'C')
        C = val4(b)*1e-6; %uF
        gkm(b) = 2*C/dt;
    elseif strcmp(type(b), 'TL')
        Zc = val4(b);
        gkm(b) = 1/Zc;
    end
end

% Build the G matrix
for b = 2:bmax
    k = from(b);
    m = to(b);
    if m == 0
        G(k,k) = G(k,k) + gkm(b);
    elseif k == 0
        G(m,m) = G(m,m) + gkm(b);
    else
        G(k,k) = G(k,k) + gkm(b);
        G(m,m) = G(m,m) + gkm(b);
        if strcmp(type(b), 'TL') == 0
            % If the branch is NOT a transmission line then calculate
            % off-diagonals
            G(k,m) = G(k,m) - gkm(b);
            G(m,k) = G(m,k) - gkm(b);
        end
    end
end
GAA = G(1:D, 1:D);
GAB = G(1:D, D+1:N);
GBA = G(D+1:N, 1:D);
GBB = G(D+1:N, D+1:N);

%% Start time as zero
n = 1;
time = 0;
t(n) = time;

THTACtl = 1; % Enable THTA in the first step of simulation

%% Simulation Loop
fprintf('-------------------- Simulation Loop ... --------------\n');

% Progress display 10%
np = fix(npoints/10);
fprintf(1, 'Progress:  00 %%\n');

while n <= npoints
       
    %% Check Switches for changes:
    % Have any switch position changed? Or is this the first step (dt)?
    % If YES, then alter matrix G for specific switch positions
    for b = 2:bmax
        if ( strcmp(type(b),'S') )
            tclose = val5(b); %tclose = val5(b) - eps;
            topen  = val6(b); %topen = val6(b)  - eps
            swchg  = 0;
            idxsw  = val8(b);
            idxplt = val9(b);
            
            if ( (SwOper(idxsw, 1) == 0) && (time >= tclose) && (val4(b) == RswOpen) )
                % This switch meets all the criterias to be closed
                
                fprintf('SWITCH %d closed at t = %2.50f.\n', b, time);
                
                val4(b) = RswClosed;
                sw_status(idx) = 1;
                SwOper(idxsw, 1) = 1;
                swchg = 1;
                sw_change_g = - 1/RswOpen + 1/RswClosed;    % add 1/RswClosed, rem 1/RswOpen
            elseif ( (SwOper(idxsw, 2) == 0) && (time >= topen) && (val4(b) == RswClosed) && ...
                    ((abs(ikm(idxplt, n-0)) <= val7(b)) || (ikm(idxplt, n-1) * ikm(idxplt, n-0) < 0))  )
                % This switch meets all the criterias to be opened
                
                fprintf('SWITCH %d opened at t = %2.50f.\n', b,time);
                
                val4(b) = RswOpen;
                sw_status(idx) = 0;
                SwOper(idxsw, 2) = 1;
                swchg = 1;
                sw_change_g = - 1/RswClosed + 1/RswOpen;    % add 1/RswOpen, rem 1/RswClosed
            end
            
            if swchg == 1
                for k = 1:N
                    node_pointer(k) = k;
                end
                % At least one switch changed its state.
                GAA = G(1:D, 1:D);
                GAB = G(1:D, D+1:N);
                GBA = G(D+1:N, 1:D);
                GBB = G(D+1:N, D+1:N);
                THTACtl = 1;
                
                R = val4(b);
                gkm(b) = 1/R;
                k = from(b);
                m = to(b);
                if m == 0
                    G(k,k) = G(k,k) + sw_change_g;
                elseif k == 0
                    G(m,m) = G(m,m) + sw_change_g;
                else
                    G(k,k) = G(k,k) + sw_change_g;
                    G(m,m) = G(m,m) + sw_change_g;
                    G(k,m) = G(k,m) - sw_change_g;
                    G(m,k) = G(m,k) - sw_change_g;
                end
            end
        end
    end
    
    %% Calculate the history sources at time t
    % to be used in t+dt for R,L,C elements and in t+Tau for TL
    Ihold = Ih;
    Ihnew = zeros(nh,1);
    for b = 2:bmax
        k = from(b);
        m = to(b);
        idxhist = val8(b);
        idxplt = val9(b);
        if strcmp(type(b), 'L')
            if k == 0
                Ihnew(idxhist) = 2 * gkm(b) * V(m,n) + Ihold(idxhist);
            elseif m == 0
                Ihnew(idxhist) = -2 * gkm(b) * V(k,n) + Ihold(idxhist);
            else
                Ihnew(idxhist) = -2 * gkm(b) * (V(k,n)-V(m,n)) + Ihold(idxhist);
            end
        elseif strcmp(type(b), 'C')
            if k == 0
                Ihnew(idxhist) = -2 * gkm(b) * V(m,n) - Ihold(idxhist);
            elseif m == 0
                Ihnew(idxhist) = 2 * gkm(b) * V(k,n) - Ihold(idxhist);
            else
                Ihnew(idxhist) = 2 * gkm(b) * (V(k,n)-V(m,n)) - Ihold(idxhist);
            end
        elseif strcmp(type(b), 'TL')
            if k == 0
                Ihkm(idxhist,n) = 2*gkm(b)*V(m,n)-IhmkTAU(idxhist);
                Ihmk(idxhist,n) = -IhkmTAU(idxhist);
            elseif m == 0
                Ihkm(idxhist,n) = -IhmkTAU(idxhist);
                Ihmk(idxhist,n) = 2*gkm(b)*V(k,n)-IhkmTAU(idxhist);
            else
                Ihkm(idxhist,n) = 2*gkm(b)*V(m,n)-IhmkTAU(idxhist);
                Ihmk(idxhist,n) = 2*gkm(b)*V(k,n)-IhkmTAU(idxhist);
            end
        end
    end
    
    %% THTA Control Here!
    if Damp == 1
        THTACtl = 0;
    end
    if THTACtl > 0
        if THTACtl == 1
            fprintf('THTA activated at t   = %2.50f.\n', time);
        end
        
        if THTACtl < 3
            fprintf('THTA step %d at t      = %2.50f.\n', THTACtl, time);
            
            Ih = (Ih + Ihnew)/2;
            
            THTACtl = THTACtl + 1;
            time = time + dt/2;
        else
            Ih = Ihnew;
            
            THTACtl = 0;
            time = time + dt;
        end
    else
        Ih = Ihnew;        
        % Regular operation
        % Increment the time
        time = time + dt;
    end

    % Store the time in the vector t
    n = n + 1;
    t(n) = time;
    
    %% Build the VB vector for the time t
    x = 1;
    for b = 2:bmax
        k = from(b);
        m = to(b);
        if strcmp(type(b), 'EDC')
            % DC Voltage Source
            VB(x,1)= val4(b);
            x = x+1;
        elseif strcmp(type(b), 'EAC')
            % AC Voltage Source
            VB(x,1)= val4(b)*cos(2*pi*val6(b)*time + val5(b)*(pi/180) );
            x = x+1;
        elseif strcmp(type(b), 'ETRI')==1
            % Triangular Voltage Source
            Fpeak = val4(b);
            tpeak = val5(b);
            t50 = val6(b);
            if time <= tpeak
                a = Fpeak/tpeak;
                Es = a*time;
            else
                a = -(Fpeak/2)/(t50 - tpeak);
                b1 = (Fpeak/2) - a*t50;
                Es = a*time + b1;
                t0 = -b1/a;
                if Es < 0
                    Es = 0;
                end
            end
            VB(x,1) = Es;
            x = x + 1;
        end
    end
    
    %% Clear vector I
    I = zeros(N,1);
    
    %% Add Current Sources evaluated at time t into vector I
    for b = 2:bmax
        k = from(b);
        m = to(b);
        Is = 0;
        if strcmp(type(b), 'IDC')
            % DC Current Source
            Is = val4(b);
        elseif strcmp(type(b), 'IAC')
            % AC Current Source
            Is = val4(b)*cos(2*pi*val6(b)*time + val5(b)*(pi/180) );
        elseif strcmp(type(b), 'ITRI')
            % Triangular Current Source
            Fpeak = val4(b);
            tpeak = val5(b);
            t50 = val6(b);
            if t(n) <= tpeak
                a = Fpeak/tpeak;
                Is = a*t(n);
            else
                a = -(Fpeak/2)/(t50 - tpeak);
                b1 = (Fpeak/2) - a*t50;
                Is = a*time + b1;
                if Is < 0
                    Is = 0;
                end
            end
        end
        if m == 0;
            I(k)= I(k) + Is;
        elseif k == 0;
            I(m)= I(m) - Is;
        else
            I(k)= I(k) + Is;
            I(m)= I(m) - Is;
        end
        Ifc(b) = Is;
    end
    
    %% Add History Current Sources, evaluated at time t-dt for R,L,C elements
    % or evaluated at t-Tau for TL, into vector I
    %
    for b = 2:bmax
        k = from(b);
        m = to(b);
        % idx = transmission line index for history terms
        % or the lumped element index for history term
        idx = val8(b);
        if strcmp(type(b), 'L') || strcmp(type(b), 'C')
            if m == 0;
                I(k)= I(k) + Ih(idx);
            elseif k == 0;
                I(m)= I(m) - Ih(idx);
            else
                I(k)= I(k) + Ih(idx);
                I(m)= I(m) - Ih(idx);
            end
        elseif strcmp(type(b), 'TL') && n - val6(b) > 1
            IhmkTAU(idx) = (Ihmk(idx,n-val6(b)-1)+val7(b)/dt*(Ihmk(idx,n-val6(b))-Ihmk(idx,n-val6(b)-1)));
            IhkmTAU(idx) = (Ihkm(idx,n-val6(b)-1)+val7(b)/dt*(Ihkm(idx,n-val6(b))-Ihkm(idx,n-val6(b)-1)));
            if k == 0
                I(m)= I(m) + IhmkTAU(idx);
            elseif m == 0
                I(k)= I(k) + IhkmTAU(idx);
            else
                I(k)= I(k) + IhkmTAU(idx);
                I(m)= I(m) + IhmkTAU(idx);
            end
        end
        
    end
    
    %% Build vector IA for the time t
    IA = I(1:D, 1);
    
    %% Build the vector RHSA for the time t
    RHSA = IA - GAB*VB;
    
    %% Solve for vector VA at the time t
    VA = GAA\RHSA;
    IB = GBA*VA+GBB*VB;
    I = [IA; IB];
    %% Build vector V at the time t (i.e.,  time counter or point n)
    V(:, n) = [VA; VB];
    
    %% Calculate required output plotting for branch voltages and branch currents at time t
    for b = 2:bmax
        if plt(b) > 0
            k = from(b);
            m = to(b);
            % identifying the branch number which have required output
            % plotting
            idx = val8(b);
            x = val9(b);
            if m == 0
                vkm(x,n) = V(k,n);
            elseif k == 0
                vkm(x,n) = - V(m,n);
            else
                vkm(x,n) = V(k,n) - V(m,n);
            end
            if strcmp(type(b), 'R')
                ikm(x,n) = gkm(b)*vkm(x,n);
                pkm(x,n) = vkm(x,n)*ikm(x,n);
                ekm(x,n) = ekm(x,n-1) + (dt/2)*( pkm(x,n-1) + pkm(x,n) );
            elseif strcmp(type(b), 'L')
                ikm(x,n) = gkm(b)*vkm(x,n) - Ih(idx);
                pkm(x,n) = vkm(x,n)*ikm(x,n);
                ekm(x,n) = ekm(x,n-1) + (dt/2)*( pkm(x,n-1) + pkm(x,n) );
            elseif strcmp(type(b), 'C')
                ikm(x,n) = gkm(b)*vkm(x,n) - Ih(idx);
                pkm(x,n) = vkm(x,n)*ikm(x,n);
                ekm(x,n) = ekm(x,n-1) + (dt/2)*( pkm(x,n-1) + pkm(x,n) );
            elseif strcmp(type(b), 'TL')
                ikm(x,n) = gkm(b)*V(k,n) - IhkmTAU(idx);
                imk(x,n) = gkm(b)*V(m,n) - IhmkTAU(idx);
                pkm(x,n) = V(k,n)*ikm(x,n);
                pmk(x,n) = V(m,n)*imk(x,n);
                ekm(x,n) = ekm(x,n-1) + (dt/2)*( pkm(x,n-1) + pkm(x,n) );
                emk(x,n) = emk(x,n-1) + (dt/2)*( pmk(x,n-1) + pmk(x,n) );
            else %Caso contrário = caso fontes (uma única de tensão)
                if k == 0
                    ikm(x,n) = I(m);
                    noeq = m; %Nó equacionado
                else
                    ikm(x,n) = I(k);
                    noeq = k;
                end
                for b2 = 2:bmax
                    if b2 ~= b
                        if strcmp(type(b2),'L') || strcmp(type(b2),'C')
                            if from(b2) == noeq
                                ikm(x,n) = ikm(x,n) - Ih(val8(b2));
                            elseif to(b2) == noeq
                                ikm(x,n) = ikm(x,n) + Ih(val8(b2));
                            end
                        end
                        if strcmp(type(b2),'TL')
                            if from(b2) == noeq
                                ikm(x,n) = ikm(x,n) - IhkmTAU(val8(b2));
                            elseif to(b2) == noeq
                                ikm(x,n) = ikm(x,n) - IhmkTAU(val8(b2));
                            end
                        end
                        if strcmp(type(b2),'IDC') || strcmp(type(b2),'IAC') || strcmp(type(b2),'ITRI')
                            if from(b2) == noeq
                                ikm(x,n) = ikm(x,n) - Ifc(b2);
                            elseif to(b2) == noeq
                                ikm(x,n) = ikm(x,n) + Ifc(b2);
                            end
                        end
                    end
                end
                if k ~= 0
                    ikm(x,n) = -ikm(x,n);
                end
            end
        end
        if strcmp(type(b), 'S')
            ikm(x,n) = gkm(b)*vkm(x,n);
        end
    end
    
    % Plot the
    if ( n >= np )
        np = np + fix(npoints/10);
        fprintf(1, 'Progress: %3d %%\n',  round(100*n/npoints));
    end
    
    
end

fprintf('---------------- Simulation loop concluded. -----------\n');
fprintf('---------------- Ready to plot the results ... --------\n');
%% Plot the results according to the user input preferences
% Plot the Nodal Voltages

if PlotNodes(1) == 0
    fprintf('Node number for nodal voltage plotting\n');
    fprintf('"99" for all; "0" to finish\r\n');
    x = 1;
    while x <= N
        node = input('Input node number: ');
        if node == 99
            PlotNodes = 1:N;
            break;
        elseif node == 0
            break;
        elseif node ~= 0
            PlotNodes(x) = node;
        end
        x = x + 1;
    end
end

if PlotNodes(1) > 0
    PlotNodes = unique(PlotNodes);
    for x = 1:length(PlotNodes)
        if PlotNodes(x) ~= 0
            figure;
            node = PlotNodes(x);
            plot(t, V(node, :), 'b-', 'LineSmoothing', 'on');
            title(['Nodal Voltage (',num2str(node), ')']);
            xlabel('time [s]');
            ylabel('voltage [V]');
            xlim(tmax*[-0.001 1]);
            grid;
        end
    end
end

for b  = 2:bmax
    k = from(b);
    m = to(b);
    p = val9(b);
    if plt(b) == 1
        figure;
        plot(t,ikm(p, :), 'b-', 'LineSmoothing', 'on');
        title( strcat( type(b), num2str(b), ' Branch Current  i(',num2str(k), ', ',num2str(m), ')' ));
        xlabel('time [s]');
        ylabel('current [A]');
        grid;
        xlim(tmax*[-0.001 1]);
        
        if strcmp(type(b), 'TL')
            hold on;
            plot(t,imk(p, :), 'r-.', 'LineSmoothing', 'on');
            title(strcat( type(b), num2str(b), ' Branch Current TL (',num2str(k), ', ',num2str(m), ')' ));
        end
        
    elseif plt(b) == 2
        figure;
        plot(t,vkm(p, :), 'b-', 'LineSmoothing', 'on');
        title( strcat( type(b), num2str(b), ' Branch Voltage  v(',num2str(k), ', ',num2str(m), ')' ));
        xlabel('time [s]');
        ylabel('voltage [V]');
        grid;
        xlim(tmax*[-0.001 1]);
        if strcmp(type(b), 'TL')
            hold on;
            plot(t,imk(p, :), 'r-.', 'LineSmoothing', 'on');
            title( strcat( type(b), num2str(b), ' Branch Current TL (',num2str(k), ', ',num2str(m), ')' ));
        end
    elseif plt(b) == 3
        figure;
        subplot(2,1,1),  plot(t,vkm(p, :), 'b-', 'LineSmoothing', 'on');
        title( strcat( type(b), num2str(b), ' Branch Voltage  v(',num2str(k), ', ',num2str(m), ')' ));
        xlabel('time [s]');
        ylabel('voltage [V]');
        grid;
        xlim(tmax*[-0.001 1]);
        subplot(2,1,2),  plot(t,ikm(p, :), 'r-', 'LineSmoothing', 'on');
        title( strcat( type(b), num2str(b), ' Branch Current  i(',num2str(k), ', ',num2str(m), ')' ));
        xlabel('time [s]');
        ylabel('current [A]');
        grid;
        xlim(tmax*[-0.001 1]);
        if strcmp(type(b), 'TL')
            hold on;
            plot(t,imk(p, :), 'r-.', 'LineSmoothing', 'on');
            title( strcat( type(b), num2str(b), ' Branch Current TL (',num2str(k), ', ',num2str(m), ')' ));
        end
    elseif plt(b) == 4
        figure;
        subplot(2,1,1),  plot(t,pkm(p, :), 'm-', 'LineSmoothing', 'on');
        title( strcat( type(b), num2str(b), ' Branch Power  p(',num2str(k), ', ',num2str(m), ')' ));
        xlabel('time [s]');
        ylabel('power [W]');
        grid;
        xlim(tmax*[-0.001 1]);
        if strcmp(type(b), 'TL')
            hold on;
            plot(t,pmk(p, :), 'm-.', 'LineSmoothing', 'on');
            title( strcat( type(b), num2str(b), ' Branch Power TL (',num2str(k), ', ',num2str(m), ')' ));
        end
        subplot(2,1,2),  plot(t,ekm(p, :), 'g-', 'LineSmoothing', 'on');
        title( strcat( type(b), num2str(b), ' Branch Energy  e(',num2str(k), ', ',num2str(m), ')' ));
        xlabel('time [s]');
        ylabel('energy [J]');
        grid;
        xlim(tmax*[-0.001 1]);
        if strcmp(type(b), 'TL')
            hold on;
            plot(t,emk(p, :), 'g-.', 'LineSmoothing', 'on');
            title( strcat( type(b), num2str(b), ' Branch Energy TL (',num2str(k), ', ',num2str(m), ')' ));
        end
    elseif plt(b) == 5
        figure;
        subplot(2,2,1),  plot(t,vkm(p, :), 'b-', 'LineSmoothing', 'on');
        title( strcat( type(b), num2str(b), ' Branch Voltage  v(',num2str(k), ', ',num2str(m), ')' ));
        xlabel('time [s]');
        ylabel('voltage  [V]');
        grid;
        xlim(tmax*[-0.001 1]);
        subplot(2,2,3),  plot(t,ikm(p, :), 'r-', 'LineSmoothing', 'on');
        title( strcat( type(b), num2str(b), ' Branch Current  i(',num2str(k), ', ',num2str(m), ')' ));
        xlabel('time [s]');
        ylabel('current [A]');
        grid;
        xlim(tmax*[-0.001 1]);
        if strcmp(type(b), 'TL')
            hold on;
            plot(t,imk(p, :), 'r-.', 'LineSmoothing', 'on');
            title( strcat( type(b), num2str(b), ' Branch Current TL (',num2str(k), ', ',num2str(m), ')' ));
        end
        subplot(2,2,2),  plot(t,pkm(p, :), 'm-', 'LineSmoothing', 'on');
        title( strcat( type(b), num2str(b), ' Branch Power  p(',num2str(k), ', ',num2str(m), ')' ));
        xlabel('time [s]');
        ylabel('power [W]');
        grid;
        xlim(tmax*[-0.001 1]);
        if strcmp(type(b), 'TL')
            hold on;
            plot(t,pmk(p, :), 'm-.', 'LineSmoothing', 'on');
            title( strcat( type(b), num2str(b), ' Branch Power TL (',num2str(k), ', ',num2str(m), ')' ));
        end
        subplot(2,2,4),  plot(t,ekm(p, :), 'g-', 'LineSmoothing', 'on');
        title( strcat( type(b), num2str(b), ' Branch Energy e(',num2str(k), ', ',num2str(m), ')' ));
        xlabel('time [s]');
        ylabel('energy [J]');
        grid;
        xlim(tmax*[-0.001 1]);
        if strcmp(type(b), 'TL')
            hold on;
            plot(t,emk(p, :), 'g-.', 'LineSmoothing', 'on');
            title( strcat( type(b), num2str(b), ' Branch Energy TL (',num2str(k), ', ',num2str(m), ')' ));
        end
    end
end

width  = 5;    % Largura, em inches
height = 5;    % Altura, em inches

width = width*100;
height = height*100;

b = 1;
fhs = flipud(findobj('Type', 'figure'));

screensize = get(0,'ScreenSize');
dsz = (screensize(3) - width)/(length(fhs)-1);

for b = length(fhs):-1:1
    fh = fhs(b);
    figure(fh);
    pos = get(fh, 'Position');
    set(gcf, 'Position', [floor((b-1)*dsz) 0 width height]);    %<- Set size
    
    lim = ylim();
    dy = (lim(2) - lim(1)) * 0.001; % Add 1% limit
    ylim([ (lim(1) - dy) (lim(2) + dy) ]);
end