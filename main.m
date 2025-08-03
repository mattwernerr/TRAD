clear, clc

datadir = "data/data/data2/";

% Get all the files within the ./data/data/data2 directory
files = dir(datadir);

% Collect all the filenames as a vector of strings
for k = size(files,1):-1:1
    allFileNames(k,1) = string(files(k).name);
end

% Iterate through each of the files
for fileNumber = 1:size(files,1)
    % Get the name of this file
    fileName = files(fileNumber).name;
    
    % Store the name of this file
    stats.filename(fileNumber,1) = string(fileName);
    
    % Skip this iteration if the file doesn't  contain ".tsv"
    if (~contains(fileName, ".tsv"))
        continue
    end
    
    % Read the data
    [t, X, q] = readData(datadir + fileName);

    % Apply fixes to the data
    X = fixX(X);
    q = fixq(q);

    %%%
    % Plot the data
    % plotData(t,X,q)
    %%%

    % Calculate the force and torque using a varying number of window
    % sizes, from 1 to 30. Note that window sizes 1 to 3 are all relatively
    % the same for the Savitzky-Golay filter (sgolay).
    for windowSize = 1:30
        % Take the first derivatives (Xdot, qdot), which are unfiltered.
        % Also take the second derivatives (Xddot, qddot), which are filtered.
        [Xdot, qdot, Xddot, qddot] = XqDerivatives(X,q,t,'sgolay',windowSize);

        % Calculate the angular velocity of the body frame relative to the inertial
        % frame, expressed in the body frame
        wb = angularVelocity(q,qdot);

        % Since wb uses qdot, it might be a little noisy so apply a filter
        wb = smoothdata(wb, 1, 'sgolay', windowSize);

        % Take a time derivative to get the angular acceleration of the body frame
        % relative to the inertial frame, expressed in the body frame AND the
        % inertial frame. This is the only vector that does this.
        wbdot = takeDerivative(wb, t);

        %% Input parameters. Ensure units are MKS
        g2kg = 1/1000;
        cm2m = 1/100;
        mm2m = 1/1000;
        in2m = 0.0254;

        % Mass
        mass.HWD = 802.114  * g2kg;
        mass.head = 5000    * g2kg;

        % Center of mass
        r.FUSION360_to_HWDforehead.CS_FUSION360 = [3.1; 6; 12]           * cm2m;
        r.FUSION360_to_HWDCOM.CS_FUSION360 = [116.041; 124.409; 72.998]  * mm2m;

        r.HWD_to_ANATOMIC.CS_body = [-(2/3)*8.5; 0; -7.25+2.75] * in2m; % GUESS!!!
        r.ANATOMIC_to_headCOM.CS_body = [0.08; -0.08; 2.77] * cm2m; % Original units: cm

        % Moment of inertia
        I.HWD.HWDCOM.CS_FUSION360 = [ 3.414 ,  0.2426, -0.7053;
                                   0.2436,  6.128 ,  0.1035;
                                  -0.7053,  0.1035,  7.670 ]*1e6 * (g2kg * mm2m^2); % Original units: g*mm2

        % The anatomic frame is aligned with the body frame
        I.head.headCOM.CS_anatomic = [109.43 0 0; 0 148.44 0; 0 0 135.88]*(5.00/3.30) * cm2m^2; % Original units: kg*cm2

        %% Transform everything to the body frame
        % Transformation matrix from FUSION360's frame to the body frame
        Rbody_F360 = [-1,0,0;0,-1,0;0,0,1];

        % Center of mass
        r.FUSION360_to_HWDforehead.CS_body = Rbody_F360*r.FUSION360_to_HWDforehead.CS_FUSION360;
        r.FUSION360_to_HWDCOM.CS_body = Rbody_F360*r.FUSION360_to_HWDCOM.CS_FUSION360;

        % Moment of inertia
        I.HWD.HWDCOM.CS_body = Rbody_F360*I.HWD.HWDCOM.CS_FUSION360*(Rbody_F360');

        %% Calculate the combined HWD+head rigid body
        % Mass
        mass.total = mass.head + mass.HWD; % kg

        % Center of mass
        r.HWDforehead_to_HWDCOM.CS_body  = r.FUSION360_to_HWDCOM.CS_body - r.FUSION360_to_HWDforehead.CS_body;
        r.HWDforehead_to_headCOM.CS_body = r.HWD_to_ANATOMIC.CS_body     + r.ANATOMIC_to_headCOM.CS_body;
        r.total.COM.CS_body = (mass.HWD*r.HWDforehead_to_HWDCOM.CS_body + mass.head*r.HWDforehead_to_headCOM.CS_body)/mass.total;

        % Moment of inertia
        I.HWD.totalCOM.CS_body = I.HWD.HWDCOM.CS_body + mass.HWD*(r.HWDforehead_to_HWDCOM.CS_body'*r.HWDforehead_to_HWDCOM.CS_body*eye(3) - r.HWDforehead_to_HWDCOM.CS_body*(r.HWDforehead_to_HWDCOM.CS_body'));
        I.head.totalCOM.CS_body = I.head.headCOM.CS_anatomic + mass.HWD*(r.HWDforehead_to_headCOM.CS_body'*r.HWDforehead_to_headCOM.CS_body*eye(3) - r.HWDforehead_to_headCOM.CS_body*(r.HWDforehead_to_headCOM.CS_body'));
        if (mass.head > 0)
            I.total.totalCOM.CS_body = I.HWD.totalCOM.CS_body + I.head.totalCOM.CS_body;
        else
            I.total.totalCOM.CS_body = I.HWD.totalCOM.CS_body;
        end

        %% Write the total values of the HWD+head system as simple variables
        % Note that everything is written in the body frame (constants!)
        %
        % The total mass
        mtot = mass.total;
        %
        % Relative position of the center of gravity, c = r_C - r_b, where r_b = 0
        % since we're writing the c.g. with respect to the body frame which by
        % definition is located at its own origin
        ctot = r.total.COM.CS_body;
        %
        % Moment of inertia
        Ibtot = I.total.totalCOM.CS_body;

        %% Calculate the force (F) and torque (M) about the HWD's sensor
        [F, M] = ForceNTorque(mtot, ctot, Ibtot, Xddot, wb, wbdot, q);

        % Keep data for calculations after the first 2 seconds but before 4:00.
        % Also remove the last trailing bits of data in case the run is under 4
        % minutes, which hopefully covers any spikes at the end. Set the data
        % values to 0 to keep numerical values while not contributing to norms.
        invalid_t_idx = find(~(2 < seconds(t) & minutes(t) < 4 & seconds(t) < seconds(t(end-5))));
        F(invalid_t_idx, :) = 0;
        M(invalid_t_idx, :) = 0;

        % Find where the time jumps
        t_invalid = t(invalid_t_idx);
        dt_invalid = seconds(diff(t_invalid));
        valid_duration = dt_invalid(dt_invalid > mean(dt_invalid));

        %% Smooth the results
        % Smooth the results
        smoothed_F = smoothdata(F, 1, 'sgolay', windowSize);
        smoothed_M = smoothdata(M, 1, 'sgolay', windowSize);

        %% Statistical calculations
        % Take the filtered time series and compute the RMS
        Fmag = vecnorm(smoothed_F, 2, 2);
        Mmag = vecnorm(smoothed_M, 2, 2);
                
        % Function RMS
        stats.F.L2(fileNumber, windowSize) = LpNorm(seconds(t), Fmag, 2)/sqrt(valid_duration);
        stats.M.L2(fileNumber, windowSize) = LpNorm(seconds(t), Mmag, 2)/sqrt(valid_duration);

        % Vector RMS
        stats.F.RMS(fileNumber, windowSize) = rms(Fmag);
        stats.M.RMS(fileNumber, windowSize) = rms(Mmag);
        
        % Standard deviation
        stats.F.STD(fileNumber, windowSize) = std(Fmag);
        stats.M.STD(fileNumber, windowSize) = std(Mmag);
        
        % Max
        stats.F.max(fileNumber, windowSize) = max(Fmag);
        stats.M.max(fileNumber, windowSize) = max(Mmag);
        
        % Compute the histogram of F and M.
        [stats.F.hist.counts{fileNumber, windowSize}, stats.F.hist.bins{fileNumber, windowSize}] = hist(Fmag, 50);
        [stats.M.hist.counts{fileNumber, windowSize}, stats.M.hist.bins{fileNumber, windowSize}] = hist(Mmag, 50);
        % % Compute the cumulative distribution function
        stats.F.cdf{fileNumber, windowSize} = cumsum(stats.F.hist.counts{fileNumber, windowSize}) / sum(stats.F.hist.counts{fileNumber, windowSize});
        stats.M.cdf{fileNumber, windowSize} = cumsum(stats.M.hist.counts{fileNumber, windowSize}) / sum(stats.M.hist.counts{fileNumber, windowSize});
    end
end

%% Sort data by indices
% Separate out by trial name
stats.tests{1,1} = "LHD";
stats.tests{1,2} = find(contains(allFileNames, "LHD"));

stats.tests{2,1} = "LHF";
stats.tests{2,2} = find(contains(allFileNames, "LHF"));

stats.tests{3,1} = "LFD";
stats.tests{3,2} = find(contains(allFileNames, "LFD"));

stats.tests{4,1} = "LFF";
stats.tests{4,2} = find(contains(allFileNames, "LFF"));

stats.tests{5,1} = "SHD";
stats.tests{5,2} = find(contains(allFileNames, "SHD"));

stats.tests{6,1} = "SHF";
stats.tests{6,2} = find(contains(allFileNames, "SHF"));

stats.tests{7,1} = "SFD";
stats.tests{7,2} = find(contains(allFileNames, "SFD"));

stats.tests{8,1} = "SFF";
stats.tests{8,2} = find(contains(allFileNames, "SFF"));

% Separate out by participant number
for k = 1:10
    stats.participant{k,1} = k;
    stats.participant{k,2} = find(contains(allFileNames, "P"+num2str(k)+"R"));
end

% Separate out by run number
for k = 1:8
    stats.run{k,1} = k;
    stats.run{k,2} = find(contains(allFileNames, "R"+num2str(k)+"_"));
end